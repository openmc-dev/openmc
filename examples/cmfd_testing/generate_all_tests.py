import numpy as np
import os
import sys

def create_timing_test(problem_size, problem_type, num_threads, root_dir, varying_input="none", run_openmc=False):
    num_inactive = str(15)
    num_batches = str(20)
    cmfd_begin = str(5)
    num_particles = str(100000)

    if problem_type == "1d1grp":
      upper_right = np.array([problem_size/2.0, 1.0, 1.0])
      bc = ["vacuum", "vacuum", "reflective", "reflective", "reflective", "reflective"]
      albedo = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    elif problem_type == "2d1grp":
      upper_right = np.array([problem_size/2.0, problem_size/2.0, 1.0])
      bc = ["vacuum", "vacuum", "vacuum", "vacuum", "reflective", "reflective"]
      albedo = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
    elif problem_type == "3d1grp":
      upper_right = np.array([problem_size/2.0, problem_size/2.0, problem_size/2.0])
      bc = ["vacuum", "vacuum", "vacuum", "vacuum", "vacuum", "vacuum"]
      albedo = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    lower_left = -1.0 * upper_right

    if varying_input == "size":
      # Dir name created based on varying size of problem
      dir1 = "size{0}-{1}-python-vec".format(problem_size, problem_type)
      dir2 = "size{0}-{1}-python-unvec".format(problem_size, problem_type)
      dir3 = "size{0}-{1}-fortran".format(problem_size, problem_type)

    elif varying_input == "threads":
      # Dir name created based on varying number of threads
      dir1 = "threads{0}-{1}-python-vec".format(num_threads, problem_type)
      dir2 = "threads{0}-{1}-python-unvec".format(num_threads, problem_type)
      dir3 = "threads{0}-{1}-fortran".format(num_threads, problem_type)

    elif varying_input == "none":
      dir1 = "{0}-python-vec".format(problem_type)
      dir2 = "{0}-python-unvec".format(problem_type)
      dir3 = "{0}-fortran".format(problem_type)

    # Create settings.xml files
    with open('base/settings_template.xml', 'r') as file:
      settings_file=file.read()
    dims = ((upper_right-lower_left)/2).astype(int)
    settings_file = settings_file.replace('{ent_mesh_dim}', "{}".format(" ".join(str(dim) for dim in dims)))
    settings_file = settings_file.replace('{ent_mesh_ll}', "{}".format(" ".join(str(ll) for ll in lower_left)))
    settings_file = settings_file.replace('{ent_mesh_ur}', "{}".format(" ".join(str(ur) for ur in upper_right)))
    settings_file = settings_file.replace('{source}', "{} {}".format(" ".join(str(ll) for ll in lower_left), " ".join(str(ur) for ur in upper_right)))
    settings_file = settings_file.replace('{num_inactive}', num_inactive)
    settings_file = settings_file.replace('{num_batches}', num_batches)
    settings_file = settings_file.replace('{num_particles}', num_particles)

    settings_file1 = settings_file.replace('{run_cmfd}', "false")
    settings_file2 = settings_file.replace('{run_cmfd}', "true")

    # Create geometry.xml_file
    with open('base/geometry_template.xml', 'r') as file:
      geometry_file=file.read()
    geometry_file = geometry_file.replace('{x_upper}', str(upper_right[0]))
    geometry_file = geometry_file.replace('{y_upper}', str(upper_right[1]))
    geometry_file = geometry_file.replace('{z_upper}', str(upper_right[2]))
    geometry_file = geometry_file.replace('{x_lower}', str(lower_left[0]))
    geometry_file = geometry_file.replace('{y_lower}', str(lower_left[1]))
    geometry_file = geometry_file.replace('{z_lower}', str(lower_left[2]))
    geometry_file = geometry_file.replace('{x_lower_bc}', str(bc[0]))
    geometry_file = geometry_file.replace('{x_upper_bc}', str(bc[1]))
    geometry_file = geometry_file.replace('{y_lower_bc}', str(bc[2]))
    geometry_file = geometry_file.replace('{y_upper_bc}', str(bc[3]))
    geometry_file = geometry_file.replace('{z_upper_bc}', str(bc[4]))
    geometry_file = geometry_file.replace('{z_lower_bc}', str(bc[5]))

    # Create cmfd.xml_file
    with open('base/cmfd_template.xml', 'r') as file:
      cmfd_file=file.read()
    cmfd_file = cmfd_file.replace('{cmfd_mesh_ll}', "{}".format(" ".join(str(ll) for ll in lower_left)))
    cmfd_file = cmfd_file.replace('{cmfd_mesh_ur}', "{}".format(" ".join(str(ur) for ur in upper_right)))
    cmfd_file = cmfd_file.replace('{cmfd_mesh_dim}', "{}".format(" ".join(str(dim) for dim in dims)))
    cmfd_file = cmfd_file.replace('{albedo}', "{}".format(" ".join(str(alb) for alb in albedo)))
    cmfd_file = cmfd_file.replace('{cmfd_begin}', cmfd_begin)

    # Create run_openmc_cmfd.py files
    with open('base/run_openmc_cmfd_template.py', 'r') as file:
      py_cmfd_file=file.read()
    py_cmfd_file = py_cmfd_file.replace('{cmfd_mesh_ll}', "{}".format(", ".join(str(ll) for ll in lower_left)))
    py_cmfd_file = py_cmfd_file.replace('{cmfd_mesh_ur}', "{}".format(", ".join(str(ur) for ur in upper_right)))
    py_cmfd_file = py_cmfd_file.replace('{cmfd_mesh_dim}', "{}".format(", ".join(str(dim) for dim in dims)))
    py_cmfd_file = py_cmfd_file.replace('{albedo}', "{}".format(", ".join(str(alb) for alb in albedo)))
    py_cmfd_file = py_cmfd_file.replace('{cmfd_begin}', cmfd_begin)
    py_cmfd_file = py_cmfd_file.replace('{omp_num_threads}', str(num_threads))

    py_cmfd_file1 = py_cmfd_file.replace('{vectorized}', "True")
    py_cmfd_file2 = py_cmfd_file.replace('{vectorized}', "False")

    # Create all subdirectories
    os.system("mkdir -p {}".format(root_dir))
    os.chdir(root_dir)

    # Create all files for vectorized python CMFD
    os.system("mkdir -p {}".format(dir1))
    os.chdir(dir1)
    os.system("cp ../../base/materials.xml ./")
    with open("settings.xml", "w") as file:
      file.write(settings_file1)
    with open("geometry.xml", "w") as file:
      file.write(geometry_file)
    with open("run_openmc_cmfd.py", "w") as file:
      file.write(py_cmfd_file1)
    if run_openmc:
      print("Running openmc on directory {}".format(dir1))
      os.system("python3 run_openmc_cmfd.py > results.dat")
    os.chdir("./..")

    # Create all files for non-vectorized python CMFD
    os.system("mkdir -p {}".format(dir2))
    os.chdir(dir2)
    os.system("cp ../../base/materials.xml ./")
    with open("settings.xml", "w") as file:
      file.write(settings_file1)
    with open("geometry.xml", "w") as file:
      file.write(geometry_file)
    with open("run_openmc_cmfd.py", "w") as file:
      file.write(py_cmfd_file2)
    if run_openmc:
      print("Running openmc on directory {}".format(dir2))
      os.system("python3 run_openmc_cmfd.py > results.dat")
    os.chdir("./..")

    # Create all files for Fortran CMFD
    os.system("mkdir -p {}".format(dir3))
    os.chdir(dir3)
    os.system("cp ../../base/materials.xml ./")
    with open("settings.xml", "w") as file:
      file.write(settings_file2)
    with open("geometry.xml", "w") as file:
      file.write(geometry_file)
    with open("cmfd.xml", "w") as file:
      file.write(cmfd_file)
    if run_openmc:
      print("Running openmc on directory {}".format(dir3))
      os.system("openmc -s {} > results.dat".format(num_threads))
    os.chdir("./..")

    os.chdir("./..")


if __name__ == "__main__":
    '''
    TODO: Add -r flag

    '''
    if len(sys.argv) > 1 and sys.argv[1] == '-r':
      run = True
    else:
      run = False

    # Generate timing test cases that vary problem sizes
    problem_sizes = [20., 40., 60., 80., 100.]
    problem_types = ["1d1grp", "2d1grp", "3d1grp"]
    num_threads = 1
    for problem_size in problem_sizes:
      for problem_type in problem_types:
        create_timing_test(problem_size, problem_type, num_threads, "timing-tests", varying_input="size", run_openmc=run)

    # Generate timing test cases that vary number of threads
    num_threads = [1, 2, 4, 8, 16]
    problem_types = ["1d1grp", "2d1grp", "3d1grp"]
    problem_size = 80.
    for threads in num_threads:
      for problem_type in problem_types:
        create_timing_test(problem_size, problem_type, threads, "timing-tests", varying_input="threads", run_openmc=run)

    # Generate test cases that don't have varying input, for replication purposes
    num_threads = 4
    problem_types = ["1d1grp", "2d1grp", "3d1grp"]
    problem_size = 20.
    for problem_type in problem_types:
      create_timing_test(problem_size, problem_type, num_threads, "replication-tests", run_openmc=run)


