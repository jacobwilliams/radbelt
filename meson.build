project('radbelt', ['c', 'fortran'],
        version: '0.0.1',
        default_options : [
          'warning_level=everything',
          'buildtype=release',]
        )

add_languages('fortran')

_args = [] # Extra arguments
_deps = [] # Dependencies
# _deps += dependency('lapack')

cc = meson.get_compiler('c')

radbeltlib = library('RadBelt',
              sources: [
                'src/radbelt_kinds_module.F90' ,
                'src/trmfun.f90' ,
                'src/shellig.f90' ,
                'src/radbelt_module.f90' ,
                'src/radbelt_c_module.f90' ,
              ],
              dependencies: _deps,
              install: true,
)
