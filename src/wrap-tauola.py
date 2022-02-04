import argparse
from dataclasses import dataclass, field
import os
import re
from typing import Union


@dataclass
class FortranEntity:
    '''Data for a fortran subroutine or function'''

    name: str
    category: str
    content: Union[str, list]
    called: bool = False


def wrap(sources, includes, outfile, version=None):
    '''Wrap TAUOLA as a fortran routine'''

    if version is None:
        version = 'undefined'

    code = []

    for source in sources:
        with open(source) as f:
            text = f.read()
        code.append(text)

    # Substitute includes
    code = os.linesep.join(code)
    functions_type = None
    if includes:
        for include in includes:
            name = os.path.basename(include)
            with open(include) as f:
                content = f.read()
            if name == 'funct_declar.inc':
                functions_type = {}
                head, tail = content.split('COMPLEX', 1)
                head, tail = tail.split('REAL', 1)
                head = re.sub(r'[&\n ]', '', head)
                functions_type['COMPLEX'] = [x.upper() for x in head.split(',')]

                head, tail = tail.split('DOUBLE PRECISION', 1)
                head = re.sub(r'[&\n ]', '', head)
                functions_type['REAL'] = [x.upper() for x in head.split(',')]

                tail = re.sub(r'[&\n ]', '', tail)
                functions_type['DOUBLE PRECISION'] = [
                    x.upper() for x in tail.split(',')
                ]

                content = ''
            content = os.linesep + content + os.linesep
            code = re.sub(
                r"\n\s*include\s*'\S*" + name + r"'\s*\n", content, code
            )

    # Map routines, functions etc.
    lines = code.split(os.linesep)
    code = None
    entities = {}
    entity, previous = None, None
    for line in lines:
        if entity is None:
            if previous:
                r = re.match(r'     [$]\s*([^(]+)', line, re.IGNORECASE)
                if r:
                    content = [previous, line]
                else:
                    content = None
                line, previous = previous, None
            else:
                r = re.match(
                    r'(?m)^\s+(?:'
                    r'SUBROUTINE|'
                    r'FUNCTION|'
                    r'REAL FUNCTION|'
                    r'DOUBLE PRECISION FUNCTION|'
                    r'COMPLEX FUNCTION|'
                    r'COMPLEX[*]16 FUNCTION)'
                    r'\s+([^(]+)',
                    line,
                    re.IGNORECASE,
                )
                if r:
                    content = [line]
                elif 'SUBROUTINE' in line:
                    previous = line

            if r:
                category = (
                    'FUNCTION' if 'FUNCTION' in line.upper() else 'SUBROUTINE'
                )
                entity = FortranEntity(r[1].strip(), category, content)
        else:
            entity.content.append(line)
            if re.match(r'(?m)^\s+END\s*$', line, re.IGNORECASE):
                entities[entity.name] = entity
                entity = None

    # Remove unused routines
    def check_call(root):
        if root.called:
            return
        else:
            root.called = True

        for line in root.content[1:]:
            if (not line) or ((line[0] != ' ') and (line[0] != '\t')):
                continue
            for entity in entities.values():
                if entity.name == root.name:
                    continue
                elif re.search(entity.name + r'[(,]', line, re.IGNORECASE):
                    check_call(entity)

    check_call(entities['INIMAS'])
    check_call(entities['INITDK'])
    check_call(entities['INIPHY'])
    check_call(entities['DEKAY'])

    for k in list(entities.keys()):
        if not entities[k].called:
            del entities[k]

    # Bind common blocks
    for entity in entities.values():
        blocks = set({})
        lines = entity.content
        for i, line in enumerate(lines):
            r = re.match(
                r'(?m)^\s+COMMON\s*/\s*([^\s]+)\s*/', line, re.IGNORECASE
            )
            if r:
                blocks.add(r.group(1))

        statblock = None
        if entity.name in ('DADMEL', 'DADMMU', 'DADMRO', 'DADMAA', 'DADMKS',
                           'DADNEW'):
            statblock = f'TAUOLA_WEIGHT_{entity.name.upper()}'
            blocks.add(statblock)

        if blocks:
            started = False
            for index, line in enumerate(lines):
                if index == 0:
                    continue
                if (
                    (len(line) < 6)
                    or ((line[0] != ' ') and (line[0] != '\t'))
                    or (line[5] != ' ')
                ):
                    continue
                r = re.match(
                    r'(?m)^\s+(?:CHARACTER|COMMON|COMPLEX|DATA|'
                    r'DOUBLE|EXTERNAL|INTEGER|PARAMETER|REAL|SAVE)',
                    line,
                    re.IGNORECASE,
                )
                if r:
                    started = True
                elif started:
                    break

            if statblock:
                lines.insert(index,
                    f'      COMMON /{statblock}/ WTMAX')
                index += 1
                if statblock != 'TAUOLA_WEIGHT_DADNEW':
                    lines.insert(index,
                        '      REAL*4 WTMAX')
                    index += 1

            lines.insert(index, '!')
            index += 1
            for block in sorted(blocks):
                cname = block.lower()
                if not cname.startswith('tauola_'):
                    cname = 'tauola_' + cname
                lines.insert(
                    index,
                    f"      BIND(C,NAME='{cname}') /{block}/",
                )
                index += 1
            lines.insert(index, '!')

    # Prune deps in routines
    for entity in entities.values():
        # Look for dependencies
        deps = []
        content = os.linesep.join(entity.content)
        for depname in entities.keys():
            if depname == entity.name:
                continue
            if re.search(depname + '\s*[(]', content, re.IGNORECASE):
                # This is a direct call
                deps.append(depname)
            elif re.search(
                r'GAUS2?[(]' + depname + r'\s*,', content, re.IGNORECASE
            ):
                # This case is a function argument to GAUS or GAUS2 integrator
                deps.append(depname)

        # Remove external declarations
        lines = content.split(os.linesep)
        for i, line in enumerate(lines):
            newline = line
            if re.match(
                r'(?m)^\s+(?:COMPLEX|DOUBLE PRECISION|'
                r'EXTERNAL|REAL)[*]?[0-9]*',
                line,
                re.IGNORECASE,
            ):
                newline += '%'
                for depname in deps:
                    if (depname == entity.name) or (
                        entities[depname].category == 'SUBROUTINE'
                    ):
                        continue
                    newline = re.sub(
                        depname + r'\s*[,%]', '', newline, re.IGNORECASE
                    )

                if newline[-1] == '%':
                    newline = newline[:-1]

                if newline != line:
                    if re.search(
                        r'(?:COMPLEX|DOUBLE PRECISION|'
                        r'EXTERNAL|REAL)[*]?[0-9 ]*$',
                        newline,
                        re.IGNORECASE,
                    ):
                        newline = ''
                    else:
                        newline = re.sub(',\s*$', '', newline)

                    lines[i] = newline

        # Set missing return type
        if (entity.category == 'FUNCTION') and functions_type:
            name = entity.name.upper()
            for k, v in functions_type.items():
                if name in v:
                    for j, l in enumerate(lines):
                        if 'IMPLICIT NONE' in l:
                            break
                    lines.insert(j + 1, f'      {k} {name}')
                    break

        entity.content = os.linesep.join(lines)

    # Redirect hard STOPs
    for entity in entities.values():
        entity.content = re.sub(
            r'(?:STOP|stop)\s*\n', 'CALL TAUOLA_STOP()\n', entity.content
        )

    # Redirect printing
    for entity in entities.values():
        lines = entity.content.split(os.linesep)

        # Map formats
        formats = {}
        check_cont, label = False, None
        for i, line in enumerate(lines):
            if line and (line[0] != ' '):
                continue
            r = re.match(r'(?m)^[\s0-9]+FORMAT', line, re.IGNORECASE)
            if r:
                label = int(line[1:6])
                formats[label] = line[r.end() + 1 : -1]
                check_cont = True
                lines[i] = '!' + line[1:]
                continue
            elif check_cont:
                if (len(line) >= 6) and (line[5] != ' '):
                    formats[label] = "''"  # Suppress long lines
                    lines[i] = '!' + line[1:]
                else:
                    check_cont, label = False, None

        # Prune formats
        def prune_format(fmt):
            '''Prune a fortran format string'''
            if fmt != "''":
                m = re.findall(r"'([^']+)'", fmt)
                if m:
                    return "'" + ''.join(m).strip() + "'"
                else:
                    return "''"
            else:
                return fmt

        for label, fmt in formats.items():
            formats[label] = prune_format(fmt)

        # Substitute printing
        check_cont = False
        for i, line in enumerate(lines):
            if line and (line[0] != ' '):
                continue
            r = re.match(r'(?m)^[\s0-9]+(?:PRINT|WRITE)', line, re.IGNORECASE)
            if r:
                prefix, suffix = line[: r.end() - 5], line[r.end() :]
                command = line[r.end() - 5 : r.end()].upper()
                if command == 'WRITE':
                    m = re.match(r"\s*[(][^,]*,'[^']+'", suffix)
                    if m:
                        args = prune_format(suffix[m.end() + 1 :])
                    else:
                        m = re.match(r"\s*[(][^,]*,([^)]+)[)]", suffix)
                        label = m.group(1)
                        if label == '*':
                            args = prune_format(suffix[m.end() + 1 :])
                        else:
                            try:
                                label = int(label)
                            except ValueError:
                                label = int(label.split('=', 1)[-1])
                            args = formats[label]
                else:
                    m = re.match(r'\s*([^,]+)', suffix)
                    label = m.group(1)
                    if label == '*':
                        args = prune_format(suffix[m.end() + 1 :])
                    else:
                        label = int(label)
                        args = formats[label]

                if args != "''":
                    args = args + '//CHAR(0)'
                    if len(args) > 57 - len(prefix):
                        args = (
                            os.linesep
                            + '     $'
                            + (len(prefix) - 4) * ' '
                            + args
                        )

                lines[i] = (
                    prefix
                    + 'CALL TAUOLA_PRINT('
                    + args
                    + ')'
                    + os.linesep
                    + '!'
                    + line[1:]
                )
                check_cont = True
                continue
            elif check_cont:
                if (len(line) >= 6) and (line[5] != ' '):
                    lines[i] = '!' + line[1:]
                else:
                    check_cont = False
        entity.content = os.linesep.join(lines)

    # Patch SIGEE declaration in SIGOLD
    if 'SIGOLD' in entities:
        entity = entities['SIGOLD']
        entity.content = re.sub(
            r'(?m)^[ ]+FUNCTION SIGOLD.*?\n',
            '''
      FUNCTION SIGOLD(Q2,JNPI)
      REAL*4 SIGEE
''',
            entity.content,
        )

    # Externalise the random engine (RANMAR)
    entity = entities['RANMAR']
    entity.content = '''
      SUBROUTINE RANMAR(RVEC,LENV)
      DIMENSION RVEC(*)
      CALL TAUOLA_RANDOM(RVEC,LENV)
      END SUBROUTINE RANMAR
'''

    # Dummy frame transform since decays are in the center of mass frame
    entities['TRALO4'] = FortranEntity(
        'TRALO4',
        'ROUTINE',
        '''
      SUBROUTINE TRALO4(KTO,P,Q,AMS)
      INTEGER KTO
      REAL*4 P(4),Q(4),AMS
      END SUBROUTINE TRALO4
''',
    )

    code = [
        f'''!     ==================================================================
!     This is a C-library compliant reformating of TAUOLA
!
!     The original FORTRAN code is available from the tauolapp website:
!     https://tauolapp.web.cern.ch/tauolapp (v1.1.8, LHC).
!
!     The main library function is the `tauola_decay` routine, defined
!     below. It wraps TAUOLA internal routines as closures. In addition,
!     common blocks used by TAUOLA are defined with explicit C binding,
!     using the `tauola_` prefix.
!     ==================================================================
      SUBROUTINE TAUOLA_DECAY(KTO,HX) BIND(C)
      INTEGER KTO,ITMP
      REAL*8 HX(4)
      COMMON /IPChT/ IVER
      INTEGER        IVER
      DATA           IVER/1/
      BIND(C,NAME='tauola_ipcht') /IPChT/
      COMMON /TAUOLA_VERSION/ VERSION
      CHARACTER(len=16)::     VERSION='{version[:15]}'//CHAR(0)
      BIND(C) /TAUOLA_VERSION/
!     ==================================================================
!     External routines needed by TAUOLA
!     ==================================================================
      INTERFACE
        SUBROUTINE TAUOLA_RANDOM(RVEC,LENV) BIND(C)
          DIMENSION RVEC(*)
        END SUBROUTINE TAUOLA_RANDOM
        SUBROUTINE FILHEP(N,STATUS,PID,MF,ML,DF,DL,P,AM,PFLAG)
     &  BIND(C,NAME='tauola_filhep')
          LOGICAL PFLAG
          INTEGER N,STATUS,PID,MF,DF,DL
          REAL P(4),AM
        END SUBROUTINE FILHEP
        SUBROUTINE TAUOLA_PRINT(S) BIND(C)
          USE ISO_C_BINDING
          IMPLICIT NONE
          CHARACTER(kind=C_CHAR), dimension(*), intent(in) :: S
        END SUBROUTINE TAUOLA_PRINT
        SUBROUTINE TAUOLA_STOP() BIND(C)
        END SUBROUTINE TAUOLA_STOP
      END INTERFACE
!     ==================================================================
!     Call the TAUOLA.DEKAY routine
!     ==================================================================
      IF(KTO.EQ.-1) THEN
        ITMP=IVER ! Backup version
        IVER=1    ! Enable new currents
        CALL RCHL_PARAMETERS(IVER)
        CALL INIMAS()
        CALL INITDK()
        CALL INIPHY(1.0D-1) ! XK0=0.1 is not used
        CALL DEKAY(KTO,HX)
        IVER=ITMP ! Restore version
      ELSE
        IF (IVER.EQ.1) THEN
          CALL RCHL_PARAMETERS(IVER)
        ENDIF
        CALL DEKAY(KTO,HX)
      ENDIF
      CONTAINS
'''
    ]

    for k, v in entities.items():
        code.append(v.content)

    code.append(
        '''
      END SUBROUTINE TAUOLA_DECAY
'''
    )

    # Prune lines and remap suppressed printing
    lines = os.linesep.join(code).split(os.linesep)
    new = []
    i = 0
    for line in lines:
        line = re.sub('^\t', 8 * ' ', line)
        line = re.sub('\t', 2 * ' ', line)
        line = re.sub(r' +$', '', line)
        line = re.sub(
            r"TAUOLA_PRINT[(]''[)]",
            f"TAUOLA_PRINT('tauola.f:{i + 1}: (suppressed)'//CHAR(0))",
            line,
        )
        if line:
            new.append(line)
            i += 1
    code.append('')
    code = os.linesep.join(new)

    # Dump the result
    with open(outfile, 'w') as f:
        f.write(code)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Wrap TAUOLA source files.')
    parser.add_argument('-s', dest='sources', nargs='+', help='source files')
    parser.add_argument(
        '-i', dest='includes', nargs='*', help='include file(s)'
    )
    parser.add_argument(
        '-w', dest='outfile', default='tauola.f', help='wrapper file'
    )
    parser.add_argument('-v', dest='version', help='source version')
    args = parser.parse_args()

    wrap(args.sources, args.includes, args.outfile, args.version)
