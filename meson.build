project('erfinv', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

compiler = meson.get_compiler('cpp')

if compiler.get_id() == 'clang' or compiler.get_id() == 'armclang'
  message('Adding parameters to relax constexpr restrictions on clang')
  add_global_arguments(['-fconstexpr-steps=200000'], language : 'cpp')
endif

exe_verification = executable('verification', 'verification/main.cpp',
  install : true,
  include_directories : '.'
  )

exe_benchmark= executable('bench', ['benchmark/main.cpp', 'benchmark/thirdparty/nanobench.cpp'],
  install : true,
  include_directories : ['.', './benchmark/thirdparty']
  )

tests = [
  ['Test', 'test/test.cpp'],
]

Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

foreach t : tests
  exe = executable(t.get(0), t.get(1),
    dependencies : [Catch2],
    include_directories : '.'
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach