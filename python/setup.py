from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os

class BuildCMakeFirst(build_ext):
    def run(self):
        # Ensure CMake build is done before python extension build
        if not os.path.isfile('lib/nikfemm/build/libnikfemm.so'):
            self.build_cmake()
        super().run()

    def build_cmake(self):
        # Build CMake project
        build_dir = os.path.abspath('../lib/nikfemm/build')
        os.makedirs(build_dir, exist_ok=True)
        debug_flag = '-DCMAKE_BUILD_TYPE=Debug' if self.debug else '-DCMAKE_BUILD_TYPE=Release'
        subprocess.check_call(['cmake', debug_flag, '..'], cwd=build_dir)
        # adaptively use the number of cores
        subprocess.check_call(['make', '-j'], cwd=build_dir)
        # copy the library to the package directory
        subprocess.check_call(['cp', '../lib/nikfemm/build/libnikfemm.so', '.'])

# Define Python extension
nikfemm_module = Extension(
    'nikfemm',
    sources=['nikfemm.cpp'],
    include_dirs=['../lib/nikfemm/include'],
    libraries=['nikfemm'],
    library_dirs=['../lib/nikfemm/build'],
    extra_link_args=['-Wl,-rpath,$ORIGIN'],
)

# link libnikfemm.so to the extension
setup(
    name='nikfemm',
    version='0.1',
    description='Python bindings for nikfemm',
    ext_modules=[nikfemm_module],
    cmdclass={'build_ext': BuildCMakeFirst},
    packages=[''],
    package_data={'': ['libnikfemm.so']},
    include_package_data=True,
    debug=True,
)