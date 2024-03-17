project('erfinv', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14', 'optimization=3'])

exe = executable('verification', 'verification/main.cpp',
  install : true,
  include_directories : '.'
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