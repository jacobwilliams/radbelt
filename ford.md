project: radbelt
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/radbelt
summary: radbelt
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         protected
         private
source: true
graph: true
externalize: true
preprocessor: gfortran -E
extra_mods: iso_fortran_env: https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}
