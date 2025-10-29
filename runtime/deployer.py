import filecmp
import os
import glob
import shutil
import fnmatch

from conan.internal.cache.home_paths import HomePaths
from conan.api.output import ConanOutput
from conan.internal.loader import load_python_file
from conan.internal.errors import conanfile_exception_formatter
from conan.errors import ConanException
from conan.internal.util.files import rmdir, mkdir


def _deploy_single(dep, conanfile, output_folder, folder_name):
    new_folder = os.path.join(output_folder, folder_name)
    rmdir(new_folder)
    symlinks = conanfile.conf.get("tools.deployer:symlinks", check_type=bool, default=True)
    try:
        shutil.copytree(dep.package_folder, new_folder, symlinks=symlinks)
        if os.path.exists(os.path.join(new_folder, "bin")):
            shutil.rmtree(os.path.join(new_folder, "bin"))
    except Exception as e:
        if "WinError 1314" in str(e):
            ConanOutput().error("full_deploy: Symlinks in Windows require admin privileges "
                                "or 'Developer mode = ON'", error_type="exception")
        raise ConanException(f"full_deploy: The copy of '{dep}' files failed: {e}.\nYou can "
                             f"use 'tools.deployer:symlinks' conf to disable symlinks")
    dep.set_deploy_folder(new_folder)

    config_file = os.path.join(output_folder, f"{dep.ref.name}-config.cmake")
    libfiles = glob.glob(os.path.join(new_folder, "lib", "*"))
    libs = " ".join(
        f"${{CMAKE_CURRENT_LIST_DIR}}/{folder_name}/lib/{os.path.basename(libfile)}"
        for libfile in libfiles
    )

    with open(config_file, "w") as dst:
        dst.write(f"set({dep.ref.name}_INCLUDE_DIRS ${{CMAKE_CURRENT_LIST_DIR}}/{folder_name}/include)\n")
        dst.write(f"set({dep.ref.name}_LIBRARIES {libs})\n")




def deploy(graph, output_folder, **kwargs):
    conanfile = graph.root.conanfile
    for dep in conanfile.dependencies.values():
        if dep.package_folder is None:
            continue
        folder_name = dep.ref.name #os.path.join("full_deploy", dep.context, dep.ref.name, str(dep.ref.version))
        build_type = dep.info.settings.get_safe("build_type")
        arch = dep.info.settings.get_safe("arch")
        #if build_type:
        #    folder_name = os.path.join(folder_name, build_type)
        #if arch:
        #    folder_name = os.path.join(folder_name, arch)
        #print("install dep:", dep)
        _deploy_single(dep, conanfile, output_folder, folder_name)
