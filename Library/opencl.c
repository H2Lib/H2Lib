/* ------------------------------------------------------------
 This is the file "opencl.c" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file opencl.c
 * @author Sven Christophersen
 * @date 2015
 */
#ifdef USE_OPENCL

/* C STD LIBRARY */
#include <unistd.h>
#include <string.h>
/* CORE 0 */
/* CORE 1 */
#include "opencl.h"
#include "clsettings.h"
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/****************************************************
 * Init the OpenCL system
 ****************************************************/

struct _ocl ocl_system;

void
get_opencl_devices(cl_device_id ** devices, cl_uint * ndevices)
{
  cl_uint   num_platform;
  cl_uint   num_devices, total_devices;
  cl_int    res;

  uint      i, j;
  char      buf[255];

  total_devices = 0;

  /****************************************************
   * Check for platforms
   ****************************************************/

  res = clGetPlatformIDs(0, 0, &num_platform);
  CL_CHECK(res)
    assert(num_platform > 0);
  ocl_system.num_platforms = num_platform;

  printf("%d platforms found:\n", num_platform);

  ocl_system.platforms =
    (cl_platform_id *) allocmem(num_platform * sizeof(cl_platform_id));

  res = clGetPlatformIDs(num_platform, ocl_system.platforms, NULL);
  CL_CHECK(res)

    for (j = 0; j < num_platform; ++j) {
    res =
      clGetPlatformInfo(ocl_system.platforms[j], CL_PLATFORM_NAME, 255, buf,
			NULL);
    CL_CHECK(res)
      printf("  platform[%d]:\t\"%s\"\n", j, buf);
    res = clGetPlatformInfo(ocl_system.platforms[j], CL_PLATFORM_VERSION, 255,
			    buf, NULL);
    CL_CHECK(res)
      printf("\t\t\"%s\"\n", buf);

    /****************************************************
     * Check for devices
     ****************************************************/

    res = clGetDeviceIDs(ocl_system.platforms[j], CL_DEVICE_TYPE_GPU, 0, 0,
			 &num_devices);
    CL_CHECK(res)
      total_devices += num_devices;

    ocl_system.devices =
      (cl_device_id *) allocmem(total_devices * sizeof(cl_device_id));

    res = clGetDeviceIDs(ocl_system.platforms[j], CL_DEVICE_TYPE_GPU,
			 num_devices, ocl_system.devices, NULL);
    CL_CHECK(res)

    /****************************************************
     * Print device names
     ****************************************************/
      for (i = 0; i < num_devices; ++i) {
      res = clGetDeviceInfo(ocl_system.devices[i], CL_DEVICE_NAME, 255,
			    (void *) buf, NULL);
      CL_CHECK(res)

	printf("    device[%d]:\t\"%s\"\n", i, buf);
    }

    freemem(ocl_system.devices);
  }

  assert(total_devices > 0);

  *devices = (cl_device_id *) allocmem(total_devices * sizeof(cl_device_id));
  ocl_system.devplatforms =
    (cl_platform_id *) allocmem(total_devices * sizeof(cl_platform_id));

  total_devices = 0;
  for (j = 0; j < num_platform; ++j) {
    res = clGetDeviceIDs(ocl_system.platforms[j], CL_DEVICE_TYPE_GPU, 0, 0,
			 &num_devices);
    CL_CHECK(res)

      for (i = 0; i < num_devices; ++i) {
      ocl_system.devplatforms[total_devices + i] = ocl_system.platforms[j];
    }

    res = clGetDeviceIDs(ocl_system.platforms[j], CL_DEVICE_TYPE_GPU,
			 num_devices, *devices + total_devices, NULL);
    CL_CHECK(res)

      total_devices += num_devices;
  }

  *ndevices = total_devices;

  ocl_system.max_package_size = (8 * 1024 * 1024);
}

void
set_opencl_devices(cl_device_id * devices, cl_uint ndevices,
		   cl_uint queues_per_device)
{
  cl_context_properties context_properties[3];
  cl_command_queue_properties queue_properties;
  cl_uint   total_devices, i, j;
  cl_int    res;
  char      buf[255];

  total_devices = 0;
  for (i = 0; i < ndevices; ++i) {
    if (devices[i] != 0) {
      devices[total_devices] = devices[i];
      ocl_system.devplatforms[total_devices] = ocl_system.devplatforms[i];

      res = clGetDeviceInfo(devices[total_devices], CL_DEVICE_NAME, 255,
			    (void *) buf, NULL);
      CL_CHECK(res)

	printf("Using device \"%s\"\n", buf);

      total_devices++;
    }
  }

  assert(queues_per_device > 0);

  ocl_system.devices = devices;
  ocl_system.num_devices = total_devices;
  ocl_system.num_contexts = total_devices;

  /****************************************************
   * Create contexts
   ****************************************************/

  ocl_system.contexts =
    (cl_context *) allocmem(total_devices * sizeof(cl_context));

  for (i = 0; i < total_devices; ++i) {

    context_properties[0] = (cl_context_properties) CL_CONTEXT_PLATFORM;
    context_properties[1] =
      (cl_context_properties) ocl_system.devplatforms[i];
    context_properties[2] = (cl_context_properties) 0;

    ocl_system.contexts[i] = clCreateContext(context_properties, 1,
					     (const cl_device_id *) devices +
					     i, NULL, NULL, &res);
    CL_CHECK(res)

  }

  /****************************************************
   * Create queues
   ****************************************************/

  ocl_system.queues =
    (cl_command_queue *) allocmem(total_devices * queues_per_device *
				  sizeof(cl_command_queue));

  queue_properties = CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
    | CL_QUEUE_PROFILING_ENABLE;

  for (i = 0; i < total_devices; ++i) {
    for (j = 0; j < queues_per_device; ++j) {
      ocl_system.queues[j + i * queues_per_device] =
	clCreateCommandQueue(ocl_system.contexts[i], ocl_system.devices[i],
			     queue_properties, &res);
      CL_CHECK(res)
    }
  }

  ocl_system.queues_per_device = queues_per_device;

  freemem(ocl_system.devplatforms);
  ocl_system.devplatforms = NULL;

}

const char *
get_error_string(cl_int error)
{

  switch (error) {
    // run-time and JIT compiler errors
  case 0:
    return "CL_SUCCESS";
  case -1:
    return "CL_DEVICE_NOT_FOUND";
  case -2:
    return "CL_DEVICE_NOT_AVAILABLE";
  case -3:
    return "CL_COMPILER_NOT_AVAILABLE";
  case -4:
    return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
  case -5:
    return "CL_OUT_OF_RESOURCES";
  case -6:
    return "CL_OUT_OF_HOST_MEMORY";
  case -7:
    return "CL_PROFILING_INFO_NOT_AVAILABLE";
  case -8:
    return "CL_MEM_COPY_OVERLAP";
  case -9:
    return "CL_IMAGE_FORMAT_MISMATCH";
  case -10:
    return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
  case -11:
    return "CL_BUILD_PROGRAM_FAILURE";
  case -12:
    return "CL_MAP_FAILURE";
  case -13:
    return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
  case -14:
    return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
  case -15:
    return "CL_COMPILE_PROGRAM_FAILURE";
  case -16:
    return "CL_LINKER_NOT_AVAILABLE";
  case -17:
    return "CL_LINK_PROGRAM_FAILURE";
  case -18:
    return "CL_DEVICE_PARTITION_FAILED";
  case -19:
    return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
  case -30:
    return "CL_INVALID_VALUE";
  case -31:
    return "CL_INVALID_DEVICE_TYPE";
  case -32:
    return "CL_INVALID_PLATFORM";
  case -33:
    return "CL_INVALID_DEVICE";
  case -34:
    return "CL_INVALID_CONTEXT";
  case -35:
    return "CL_INVALID_QUEUE_PROPERTIES";
  case -36:
    return "CL_INVALID_COMMAND_QUEUE";
  case -37:
    return "CL_INVALID_HOST_PTR";
  case -38:
    return "CL_INVALID_MEM_OBJECT";
  case -39:
    return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
  case -40:
    return "CL_INVALID_IMAGE_SIZE";
  case -41:
    return "CL_INVALID_SAMPLER";
  case -42:
    return "CL_INVALID_BINARY";
  case -43:
    return "CL_INVALID_BUILD_OPTIONS";
  case -44:
    return "CL_INVALID_PROGRAM";
  case -45:
    return "CL_INVALID_PROGRAM_EXECUTABLE";
  case -46:
    return "CL_INVALID_KERNEL_NAME";
  case -47:
    return "CL_INVALID_KERNEL_DEFINITION";
  case -48:
    return "CL_INVALID_KERNEL";
  case -49:
    return "CL_INVALID_ARG_INDEX";
  case -50:
    return "CL_INVALID_ARG_VALUE";
  case -51:
    return "CL_INVALID_ARG_SIZE";
  case -52:
    return "CL_INVALID_KERNEL_ARGS";
  case -53:
    return "CL_INVALID_WORK_DIMENSION";
  case -54:
    return "CL_INVALID_WORK_GROUP_SIZE";
  case -55:
    return "CL_INVALID_WORK_ITEM_SIZE";
  case -56:
    return "CL_INVALID_GLOBAL_OFFSET";
  case -57:
    return "CL_INVALID_EVENT_WAIT_LIST";
  case -58:
    return "CL_INVALID_EVENT";
  case -59:
    return "CL_INVALID_OPERATION";
  case -60:
    return "CL_INVALID_GL_OBJECT";
  case -61:
    return "CL_INVALID_BUFFER_SIZE";
  case -62:
    return "CL_INVALID_MIP_LEVEL";
  case -63:
    return "CL_INVALID_GLOBAL_WORK_SIZE";
  case -64:
    return "CL_INVALID_PROPERTY";
  case -65:
    return "CL_INVALID_IMAGE_DESCRIPTOR";
  case -66:
    return "CL_INVALID_COMPILER_OPTIONS";
  case -67:
    return "CL_INVALID_LINKER_OPTIONS";
  case -68:
    return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
  case -1000:
    return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
  case -1001:
    return "CL_PLATFORM_NOT_FOUND_KHR";
  case -1002:
    return "CL_INVALID_D3D10_DEVICE_KHR";
  case -1003:
    return "CL_INVALID_D3D10_RESOURCE_KHR";
  case -1004:
    return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
  case -1005:
    return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
  default:
    return "Unknown OpenCL error";
  }
}

void
setup_kernels(const char *src_str, const uint num_kernels,
	      const char **kernel_names, cl_kernel ** kernels)
{
  cl_uint   num_devices = ocl_system.num_devices;
  cl_uint   num_queues = ocl_system.queues_per_device;
  const char **srcs;
  size_t    size_src_str[2];
  cl_int    res;
  cl_int    i, j, k;
  cl_program program;

  size_src_str[0] = strlen(clsettings_src);
  size_src_str[1] = strlen(src_str);
  srcs = (const char **) allocmem(2 * sizeof(char *));
  srcs[0] = clsettings_src;
  srcs[1] = src_str;

  /****************************************************
   * Create OpenCL program object.
   ****************************************************/

  *kernels =
    (cl_kernel *) allocmem(num_queues * num_kernels * num_devices *
			   sizeof(cl_kernel));

  for (k = 0; k < num_devices; ++k) {

    program = clCreateProgramWithSource(ocl_system.contexts[k], 2, srcs,
					size_src_str, &res);
    CL_CHECK(res)

    /****************************************************
     * Compile OpenCL program object.
     ****************************************************/
      res = clBuildProgram(program, 1, ocl_system.devices + k,
			   "-Werror -cl-mad-enable -cl-fast-relaxed-math"
#ifdef USE_COMPLEX
			   " -DUSE_COMPLEX"
#endif
#ifdef USE_FLOAT
			   " -DUSE_FLOAT"
#endif
			   , NULL, NULL);

    if (res != CL_SUCCESS) {
      cl_device_id dev_id = ocl_system.devices[k];
      size_t    len;
      char      buffer[204800];
      cl_build_status bldstatus;
      printf("\nError %d: Failed to build program executable [ %s ]\n", res,
	     get_error_string(res));
      res = clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_STATUS,
				  sizeof(bldstatus), (void *) &bldstatus,
				  &len);
      if (res != CL_SUCCESS) {
	printf("Build Status error %d: %s\n", res, get_error_string(res));
	exit(1);
      }
      if (bldstatus == CL_BUILD_SUCCESS)
	printf("Build Status: CL_BUILD_SUCCESS\n");
      if (bldstatus == CL_BUILD_NONE)
	printf("Build Status: CL_BUILD_NONE\n");
      if (bldstatus == CL_BUILD_ERROR)
	printf("Build Status: CL_BUILD_ERROR\n");
      if (bldstatus == CL_BUILD_IN_PROGRESS)
	printf("Build Status: CL_BUILD_IN_PROGRESS\n");
      res = clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_OPTIONS,
				  sizeof(buffer), buffer, &len);
      if (res != CL_SUCCESS) {
	printf("Build Options error %d: %s\n", res, get_error_string(res));
	exit(1);
      }
      printf("Build Options: %s\n", buffer);
      res = clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_LOG,
				  sizeof(buffer), buffer, &len);
      if (res != CL_SUCCESS) {
	printf("Build Log error %d: %s\n", res, get_error_string(res));
	exit(1);
      }
      printf("Build Log:\n%s\n", buffer);
      abort();
    }

    /****************************************************
     * Build kernels
     ****************************************************/

    for (i = 0; i < num_kernels; ++i) {
      for (j = 0; j < num_queues; ++j) {
	(*kernels)[j + k * num_queues + i * num_devices * num_queues] =
	  clCreateKernel(program, kernel_names[i], &res);
	CL_CHECK(res);
      }
    }

    clReleaseProgram(program);
  }

  freemem(srcs);
}

/****************************************************
 * Tasks and taskgroups
 ****************************************************/

ptask
new_ocltask(void *data, ptask next)
{
  ptask     t;

  assert(data != NULL);

  t = (ptask) allocmem(sizeof(task));
  t->data = data;
  t->next = next;

  return t;
}

void
del_ocltask(ptaskgroup tg, ptask t)
{
  assert(t != NULL);

  tg->cleanup_task(t->data);

  freemem(t);
}

void
del_recursive_ocltasks(ptaskgroup tg, ptask t)
{
  ptask     tnext;

  while (t != NULL) {
    tnext = t->next;
    del_ocltask(tg, t);
    t = tnext;
  }
}

ptaskgroup fill_taskgroups = NULL;
omp_lock_t fill_lock;

ptaskgroup ready_taskgroups = NULL;
omp_lock_t ready_lock;

typedef enum {
  FILL_LIST,			// List of taskgroups, that are currently being filled up with tasks.
  READY_LIST			// List of taskgroups, that are ready for execution.
} taskgrouplist;

static ptaskgroup *
get_list(taskgrouplist tgl)
{
  ptaskgroup *tg = NULL;

  switch (tgl) {
  case FILL_LIST:
    tg = &fill_taskgroups;
    break;
  case READY_LIST:
    tg = &ready_taskgroups;
    break;
  default:
    fprintf(stderr, "Unexpected list!\n");
    abort();
    break;
  }

  return tg;
}

static omp_lock_t *
get_lock(taskgrouplist tgl)
{
  omp_lock_t *lock = NULL;

  switch (tgl) {
  case FILL_LIST:
    lock = &fill_lock;
    break;
  case READY_LIST:
    lock = &ready_lock;
    break;
  default:
    fprintf(stderr, "Unexpected list!\n");
    abort();
    break;
  }

  return lock;
}

static void
register_tolist_taskgroup(ptaskgroup tg, taskgrouplist tgl)
{
  omp_lock_t *lock = get_lock(tgl);
  ptaskgroup *list, tg2;

  omp_set_lock(lock);

  list = get_list(tgl);

  assert(list != NULL);

  if (*list == NULL) {

    /****************************************************
     * If 'list' is empty, then 'tg' is set as head of the
     * list.
     ****************************************************/

    *list = tg;
    tg->next = NULL;
  }
  else {

    /****************************************************
     * Check for duplicate entries in the list.
     ****************************************************/

    tg2 = *list;
    while (tg2 != NULL) {
      assert(tg != tg2);
      tg2 = tg2->next;
    }

    /****************************************************
     * If 'list' is not empty, then 'tg' is set as new
     * head of the list and 'list' is appended to 'tg'.
     ****************************************************/

    tg->next = *list;
    *list = tg;
  }

  omp_unset_lock(lock);
}

static void
unregister_fromlist_taskgroup(ptaskgroup tg, taskgrouplist tgl)
{
  omp_lock_t *lock = get_lock(tgl);
  ptaskgroup *list;
  ptaskgroup tg2;

  omp_set_lock(lock);

  list = get_list(tgl);

  assert(tg != NULL);

  tg2 = *list;

  /****************************************************
   * If list is empty, then we can't unregister anything.
   ****************************************************/

  if (tg2 == NULL) {
    omp_unset_lock(lock);
    return;
  }

  /****************************************************
   * If Element if the head of the list, then remove it
   * and set 2nd element as new head of the list.
   ****************************************************/

  if (tg2 != NULL && tg == tg2) {
    *list = tg2->next;
    omp_unset_lock(lock);
    return;
  }

  /****************************************************
   * Traverse the list until we either reach the end
   * of the list or we find the desired element.
   ****************************************************/

  while (tg2 != NULL && tg2->next != NULL && tg2->next != tg) {
    assert(tg2->next != tg2);
    tg2 = tg2->next;
  }
  assert(tg2 != NULL);

  /****************************************************
   * If we reached the end of the list we are done.
   ****************************************************/

  if (tg2->next == NULL) {
    assert(tg != tg2);
    omp_unset_lock(lock);
    return;
  }

  /****************************************************
   * Otherwise we need to extract the element somewhere
   * in the middle of the list.
   ****************************************************/

  assert(tg2->next != NULL);
  if (tg2->next == tg) {
    tg2->next = tg2->next->next;
    omp_unset_lock(lock);
    return;
  }

  omp_unset_lock(lock);
  fprintf(stderr, "Unexpected case in 'unregister_fromlist_taskgroup'!\n");
  abort();
}

ptaskgroup
new_ocltaskgroup(task_affinity affinity, ptaskgroup next,
		 void (*merge_tasks) (ptaskgroup tg, void **data),
		 void (*cleanup_merge) (void *data),
		 void (*distribute_results) (ptaskgroup tg, void *data),
		 void (*close_taskgroup) (ptaskgroup tg),
		 void (*callback_cpu) (void *data),
		 void (*callback_gpu) (void *data),
		 size_t(*getsize_task) (void *data),
		 void (*split_task) (void *data, void ***split, uint * n),
		 void (*cleanup_task) (void *data), void *data)
{
  ptaskgroup tg;

  if (affinity == CPU_ONLY || affinity == CPU_FIRST) {
    assert(callback_cpu != NULL);
  }

  if (affinity == GPU_ONLY || affinity == GPU_FIRST) {
    assert(callback_gpu != NULL && merge_tasks != NULL
	   && distribute_results != NULL);
    assert(cleanup_merge != NULL);
  }

  assert(getsize_task != NULL);
  assert(cleanup_task != NULL);
  assert(close_taskgroup != NULL);

  tg = (ptaskgroup) allocmem(sizeof(taskgroup));
  tg->tasks = NULL;
  tg->size = 0;
  tg->getsize_task = getsize_task;
  tg->affinity = affinity;
  tg->cleanup_task = cleanup_task;
  tg->split_task = split_task;
  tg->cleanup_merge = cleanup_merge;
  tg->close_taskgroup = close_taskgroup;
  tg->callback_cpu = callback_cpu;
  tg->callback_gpu = callback_gpu;
  tg->merge_tasks = merge_tasks;
  tg->distribute_results = distribute_results;
  tg->next = next;
  tg->data = data;

  register_tolist_taskgroup(tg, FILL_LIST);

  return tg;
}

void
del_taskgroup(ptaskgroup tg)
{
  del_recursive_ocltasks(tg, tg->tasks);

  freemem(tg);
}

ptaskgroup
dequeue_taskgroup(taskgrouplist tgl)
{
  omp_lock_t *lock = get_lock(tgl);
  ptaskgroup *tg, tg2;

  omp_set_lock(lock);

  tg = get_list(tgl);
  tg2 = *tg;

  if (tg2 != NULL) {
    *tg = tg2->next;
  }

  omp_unset_lock(lock);

  return tg2;
}

static    bool
add_task(ptaskgroup tg, void *data)
{
  ptask     t;
  size_t    tsize;

  assert(tg != NULL);
  assert(data != NULL);

  /****************************************************
   * Get size of the current task and check if it fits
   * into the current taskgroup or not.
   * Return 'false' in the first and 'true' in the
   * latter case.
   ****************************************************/

  tsize = tg->getsize_task(data);

  if (tg->size + tsize < ocl_system.max_package_size) {
    t = new_ocltask(data, tg->tasks);
    tg->tasks = t;
    tg->size += tsize;
    return false;
  }
  else {
    return true;
  }
}

void
add_task_taskgroup(ptaskgroup * tg, void *data)
{
  void    **split;
  uint      n, i;

  assert(tg != NULL);
  assert(*tg != NULL);
  assert(data != NULL);

  if (add_task(*tg, data)) {
    if ((*tg)->size > 0) {
      add_taskgroup_to_scheduler(*tg);
      *tg = new_ocltaskgroup((*tg)->affinity, NULL, (*tg)->merge_tasks,
			     (*tg)->cleanup_merge, (*tg)->distribute_results,
			     (*tg)->close_taskgroup, (*tg)->callback_cpu,
			     (*tg)->callback_gpu, (*tg)->getsize_task,
			     (*tg)->split_task, (*tg)->cleanup_task,
			     (*tg)->data);
    }
    if (add_task(*tg, data)) {
      (*tg)->split_task(data, &split, &n);
      for (i = 0; i < n; ++i) {
	assert(split[i] != NULL);
	add_task_taskgroup(tg, split[i]);
      }

      (*tg)->cleanup_task(data);
      freemem(split);
    }
  }
}

static void
execute_taskgroup_cpu(ptaskgroup tg)
{
  ptask     t;

  assert(tg != NULL);

  /****************************************************
   * Execute the tasklist one by one on the CPU
   ****************************************************/

  t = tg->tasks;
  while (t != NULL) {
    tg->callback_cpu(t->data);
    t = t->next;
  }
}

static void
execute_taskgroup_gpu(ptaskgroup tg)
{
  void     *merge_data;

  assert(tg != NULL);
  assert(tg->tasks != NULL);

  /****************************************************
   * Merge tasklist to a single task
   ****************************************************/

  tg->merge_tasks(tg, &merge_data);

  /****************************************************
   * Execute task on the GPU
   ****************************************************/

  tg->callback_gpu(merge_data);

  /****************************************************
   * Distribute data to the particular tasks
   ****************************************************/

  tg->distribute_results(tg, merge_data);

  /****************************************************
   * Cleanup 'merge_data'
   ****************************************************/

  tg->cleanup_merge(merge_data);
}

/****************************************************
 * OpenCL-OpenMP scheduler
 ****************************************************/

void
add_taskgroup_to_scheduler(ptaskgroup tg)
{

  /****************************************************
   * Remove taskgroup from FILL_LIST and insert into
   * READY_LIST. Executing of the taskgroup is part of
   * the scheduler itself, who is served by the
   * READY_LIST.
   ****************************************************/

  unregister_fromlist_taskgroup(tg, FILL_LIST);
  register_tolist_taskgroup(tg, READY_LIST);
}

uint      scheduler_active = 0;

static void
flush_active_taskgroup()
{
  ptaskgroup tg, tg2, *tgp;
  omp_lock_t *lock;

  /****************************************************
   * Remove all taskgroups from FILL_LIST and insert
   * into READY_LIST. Executing of the taskgroups is
   * part of the scheduler itself, who is served by the
   * READY_LIST.
   ****************************************************/

  lock = get_lock(FILL_LIST);
  omp_set_lock(lock);

  tg = *get_list(FILL_LIST);
  while (tg != NULL) {
    tg2 = tg->next;
    tg->next = NULL;
    tg->close_taskgroup(tg);
    register_tolist_taskgroup(tg, READY_LIST);
    tg = tg2;
  }
  tgp = get_list(FILL_LIST);
  *tgp = NULL;

  omp_unset_lock(lock);
}

void
start_scheduler(uint cpu_threads, uint gpu_threads)
{
  uint      nthreads = cpu_threads + gpu_threads;
  uint      active = 0;

  scheduler_active = 1;

  if (omp_get_thread_num() == 0) {
    omp_init_lock(get_lock(FILL_LIST));
    omp_init_lock(get_lock(READY_LIST));
  }

#pragma omp barrier

  if (omp_get_thread_num() > 0) {
#pragma omp parallel num_threads(nthreads)
    {
      ptaskgroup tg;

      while ((scheduler_active == 1) || (*get_list(FILL_LIST) != NULL)
	     || (*get_list(READY_LIST) != NULL)) {

#pragma omp single
	{
	  active = 0;
	  if (scheduler_active == 0) {
	    flush_active_taskgroup();
	  }
	}

	while ((tg = dequeue_taskgroup(READY_LIST)) != NULL
	       || scheduler_active == 1 || active > 0) {
	  if (tg != NULL) {

#pragma omp atomic
	    active++;

	    switch (tg->affinity) {
	    case CPU_FIRST:
	    case CPU_ONLY:
	      execute_taskgroup_cpu(tg);
	      break;
	    case GPU_FIRST:
	      if (omp_get_thread_num() < gpu_threads) {
		execute_taskgroup_gpu(tg);
	      }
	      else {
		execute_taskgroup_cpu(tg);
	      }
	      break;
	    case GPU_ONLY:
	      if (omp_get_thread_num() < gpu_threads) {
		execute_taskgroup_gpu(tg);
	      }
	      else {
		fprintf(stderr, "HOUSTON! We have a problem!\n");
	      }
	      break;
	    default:
	      fprintf(stderr, "Unknown task affinity!\n");
	      abort();
	      break;
	    }

	    del_taskgroup(tg);

#pragma omp atomic
	    active--;
	  }
	  else {
	    if (active == 0 && scheduler_active != 1) {
	      break;
	    }

	    usleep(100);
	  }
	}

#pragma omp barrier
      }
    }

    assert(fill_taskgroups == NULL);
    assert(ready_taskgroups == NULL);
    assert(active == 0);
  }
}

void
stop_scheduler()
{
  scheduler_active = 0;

#pragma omp barrier

  if (omp_get_thread_num() == 0) {
    omp_destroy_lock(get_lock(FILL_LIST));
    omp_destroy_lock(get_lock(READY_LIST));
  }
}

#endif
