/* ------------------------------------------------------------
 This is the file "opencl.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file opencl.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef OPENCL_H_
#define OPENCL_H_

#ifdef USE_OPENCL

/** @defgroup opencl opencl
 *  @brief Management routines for OpenCL computations via task scheduler.
 *
 * @{
 */

/* C STD LIBRARY */
#include <stdio.h>
#include <assert.h>
#include <CL/cl.h>
#include <omp.h>
/* CORE 0 */
#include "basic.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/**
 * @brief Structure that contains basic OpenCL objects for arbitrary
 * computations.
 *
 * When this object is set up OpenCL platforms, devices,
 * contexts and queues, that will be used for subsequent computations, are
 * determined.
 */
struct _ocl {
  /** @brief Number of OpenCL platforms that will be used for computations. */
  cl_uint num_platforms;
  /** @brief This array contains all OpenCL platform-ids that will be used.*/
  cl_platform_id *platforms;
  /** @brief Just a temporary array of OpenCL platform-ids that is used in
   between the calls of @ref get_opencl_devices and @ref set_opencl_devices.
   */
  cl_platform_id *devplatforms;
  /** @brief Number of OpenCL devices the will be used. */
  cl_uint num_devices;
  /** @brief Array of OpenCl devices that will be used. */
  cl_device_id *devices;
  /** @brief Number of OpenCL contexts that will be used. */
  cl_uint num_contexts;
  /** @brief Array of OpenCL contexts that will be used. */
  cl_context *contexts;
  /** @brief This values determines the number of queues that will be created
   per device. */
  cl_uint queues_per_device;
  /** @brief Array of all queues for all devices that will be used. */
  cl_command_queue *queues;
  /** Determines the maximal size of each processing chunk. */
  size_t max_package_size;
};

/**
 * @brief Global variable that contains basic OpenCL objects for arbitrary
 * computations.
 */
extern struct _ocl ocl_system;

/**
 * Simple macro, that checks OpenCL error codes and prints them to @c stderr.
 */
#define CL_CHECK(res) \
    {if (res != CL_SUCCESS) {fprintf(stderr,"Error \"%s\" (%d) in file %s on line %d\n", \
        get_error_string(res), res, __FILE__,__LINE__); abort();}}

/**
 * @brief Retrieve an array of available OpenCL devices.
 *
 * @param devices Pointer to an array where the results should be stored.
 * @param ndevices Number of retrieved devices will be stored in this variable.
 */
HEADER_PREFIX void get_opencl_devices(cl_device_id **devices, cl_uint *ndevices);

/**
 * @brief Set an array of OpenCL devices to be used for computations.
 *
 * @param devices An array of OpenCL device-ids. The length should be the same
 *        as the result from @ref get_opencl_devices array. Devices that should
 *        not be used have to be set to @c 0.
 * @param ndevices Number of retrieved devices will be stored in this variable.
 * @param queues_per_device Determines the number of CPU threads that should
 *        employ a single GPU.
 */
HEADER_PREFIX void set_opencl_devices(cl_device_id *devices, cl_uint ndevices,
    cl_uint queues_per_device);

/**
 * @brief Returns a string corresponding to the error code \p error.
 *
 * @param error An error code return by some OpenCL runtime function.
 *
 * @return The Name of the OpenCL macro defining the error code.
 */
HEADER_PREFIX const char *get_error_string(cl_int error);

/**
 * @brief Reads a file specified by @p filename and compiles all OpenCL kernels
 * given by the array @p kernel_names into OpenCL kernels @p kernels.
 *
 * @param src_str Source code string.
 * @param num_kernels Number of kernels that should be compiled given by @p
 *        kernel_names.
 * @param kernel_names Array of function names that should be compiled as OpenCL
 *        kernels.
 * @param kernels Resulting array of OpenCL kernels.
 */
HEADER_PREFIX void setup_kernels(const char *src_str, const uint num_kernels,
    const char **kernel_names, cl_kernel **kernels);

/**
 * @brief This enum specifies the affinity of exceution of a @ref _taskgroup
 *        "taskgroup".
 */
typedef enum {
  CPU_ONLY, /**< CPU_ONLY A @ref _taskgroup "taskgroup" should @b only be
   executed on the cpu. */
  GPU_ONLY, /**< GPU_ONLY A @ref _taskgroup "taskgroup" should @b only be
   executed on the gpu. */
  CPU_FIRST,/**< CPU_FIRST A @ref _taskgroup "taskgroup" <b>is preferred</b> to
   be executed on the cpu, but execution on the gpu might be possible. */
  GPU_FIRST /**< GPU_FIRST A @ref _taskgroup "taskgroup" <b>is preferred</b> to
   be executed on the gpu, but execution on the cpu might be possible. */
}task_affinity;

/**
 * @brief Abbreviation for a @ref _task "struct _task".
 */
typedef struct _task task;

/**
 * @brief Abbreviation for a pointer to a @ref _task "struct _task".
 */
typedef task *ptask;

/**
 * @brief Simple representation of a task.
 *
 * The data of the a task is stored in the field @ref _task.data.
 * Because this struct is usually used in the same contexts with a
 * @ref _taskgroup "taskgroup", the code to be executed should be clear.
 * Also it contains pointer @ref _task.next to the next task in the list.
 */
struct _task {
  /** Input / output data for the current task. */
  void *data;
  /** Pointer to the next task in the list. */
  ptask next;
};

/**
 * @brief Abbreviation for a @ref _taskgroup "struct _taskgroup".
 */
typedef struct _taskgroup taskgroup;

/**
 * @brief Abbreviation for a pointer to a @ref _taskgroup "struct _taskgroup".
 */
typedef taskgroup *ptaskgroup;

/**
 * @brief A collection of tasks for the same callback function.
 *
 * Depending on the taskgroup affinity @ref _taskgroup.affinity the computation
 * can either be done on the CPU by the callback
 * @ref _taskgroup.callback_cpu or optionally on the GPU by
 * the second callback function @ref _taskgroup.callback_gpu .
 * In the latter case also the callbacks @ref _taskgroup.merge_tasks,
 * @ref _taskgroup.distribute_results and @ref _taskgroup.cleanup_merge have
 * to be specified.
 * The List of tasks to be processed is given by @ref _taskgroup.tasks and
 * the field @ref _taskgroup.size determines the current matrix elements
 * to be processed by the task list.
 * The amount of storage used by the results can be calculated by the callback
 * @ref _taskgroup.getsize_task.
 * When the scheduler calls a flush on all available @ref _taskgroup "taskgroups"
 * the callback @ref _taskgroup.close_taskgroup will be called. It is mandatory
 * to supply this callback because some @ref _taskgroup "taskgroups" might need
 * a chance to close their collection of work before enqueueing them to the
 * scheduler.
 * Also inevitable is the callback @ref _taskgroup.split_task that might split
 * tasks that are too big into smaller chunks, that again might fit into a
 * pending @ref _taskgroup.
 * Finally the callback @ref _taskgroup.cleanup_task takes care of cleaning up
 * temporary data for a given @ref _task.
 */
struct _taskgroup {
  /** The List of tasks to be processed. */
  ptask tasks;
  /** Current size of the corresponding matrix elements to be processed. */
  size_t size;
  /** Determines the affinity for computation hardware. */
  task_affinity affinity;
  /** Pointer to the next taskgroup in the global list */
  ptaskgroup next;
  /** Additional information about a taskgroup.*/
  void *data;

  /** Callback for the collection of input data from tasks for the GPU
   computation. */
  void (*merge_tasks)(ptaskgroup t, void **data);
  /** Cleanup callback for the merged data. */
  void (*cleanup_merge)(void *data);
  /** Callback for the distribution of output data from GPU computation. */
  void (*distribute_results)(ptaskgroup t, void *data);
  /** Close pending work for a taskgroup. */
  void (*close_taskgroup)(ptaskgroup tg);

  /** Callback for the computation on the CPU. */
  void (*callback_cpu)(void *data);
  /** Callback for the computation on the GPU. */
  void (*callback_gpu)(void *data);
  /** Calculate the storage amount for the results. */
  size_t (*getsize_task)(void *data);
  /** Splits a @ref _task into number of smaller tasks. The count is returned
   * by @p n, the array of new tasks is returned by @p splits.*/
  void (*split_task)(void *data, void ***splits, uint *n);
  /** Cleanup callback for each task of the taskgroup. */
  void (*cleanup_task)(void *data);
};

/**
 * @brief Create a new task.
 *
 * A new @ref _task object will be created containing the input / output data
 * of the task specified by @p data and a pointer to the next element in
 * the task list given by @p next.
 *
 * @param data Data to processed by the new task.
 * @param next Pointer to the next task in the list.
 *
 * @return The newly created task object.
 */
HEADER_PREFIX ptask new_ocltask(void *data, ptask next);

/**
 * @brief Delete a task.
 *
 * @param tg @ref _taskgroup "Taskgroup" object that contains the correct cleanup
 *        callback for the deletion of @p t.
 * @param t @ref _task "Task" to be deleted.
 */
HEADER_PREFIX void del_ocltask(ptaskgroup tg, ptask t);

/**
 * @brief Create a new list of tasks with the same code to execute.
 *
 * Create a new @ref _taskgroup object and set the callback functions in a way
 * that the correct piece of code will be executed on the data provided by the
 * particular tasks.
 *
 * @param affinity
 * @param next
 * @param merge_tasks Callback function that merges the data of the task list
 *        into a single data element for processing on the GPU.
 * @param cleanup_merge Callback function which is called after completion of
 *        @p merge_tasks to cleanup its intermediate results.
 * @param distribute_results Callback function that distributes the results of
 *        @p callback_gpu corresponding to the task list.
 * @param close_taskgroup This callback function is responsible to safely close
 *        a non-empty but not completely filled taskgroup.
 * @param getsize_task Calculate the size of the results.
 * @param callback_cpu Code to be executed on the CPU for the list of tasks.
 * @param callback_gpu Code to be executed on the GPU for the list of tasks.
 *        When @p callback_gpu is set, the callbacks @p merge_tasks and
 *        @p distribute_results have to be set correctly as well.
 * @param getsize_task The result of this callback function is the size in Byte,
 *        that a single task, belonging to the current taskgroup, has.
 * @param split_task This callback function will split a single task into a
 *        number of smaller tasks, that might fit better into a taskgroup than
 *        the original one.
 * @param cleanup_task Callback to free memory used by a task.
 * @param data Additional information for a taskgroup
 * @return Returns the newly created @ref taskgroup object.
 */
HEADER_PREFIX ptaskgroup new_ocltaskgroup(task_affinity affinity, ptaskgroup next,
    void (*merge_tasks)(ptaskgroup tg, void **data),
    void (*cleanup_merge)(void *data),
    void (*distribute_results)(ptaskgroup tg, void *data),
    void (*close_taskgroup)(ptaskgroup tg), void (*callback_cpu)(void *data),
    void (*callback_gpu)(void *data), size_t (*getsize_task)(void *data),
    void (*split_task)(void *data, void ***split, uint *n),
    void (*cleanup_task)(void *data), void *data);

/**
 * @brief Delete a @ref taskgroup object.
 *
 * @param tg @ref taskgroup object to be deleted.
 */
HEADER_PREFIX void del_taskgroup(ptaskgroup tg);

/**
 * @brief Adds a new @ref _task to a @ref _taskgroup.
 *
 * The task is enqueued into the task list of the @ref _taskgroup object @p tg.
 *
 * @param tg The @ref _taskgroup object where the task has to be enqueued to.
 * @param data The data to be enqueued.
 *
 * @attention If the function returns @a true the task is @b not added to the
 *            queue.
 */
HEADER_PREFIX void add_task_taskgroup(ptaskgroup *tg, void *data);

/**
 * @brief Enqueues a full taskgroup to the execution queue of the scheduler.
 *
 * @param tg @ref _taskgroup to be enqueued.
 */
HEADER_PREFIX void add_taskgroup_to_scheduler(ptaskgroup tg);

/**
 * @brief Starts the task scheduler with @p cpu_threads Threads on the CPU
 *        and @p gpu_threads on the GPU.
 * @param cpu_threads Number of CPU threads for task executions
 * @param gpu_threads Number of CPU threads, that employ the available GPUs
 *        for task executions.
 *
 * @attention Both, the master and the slave thread must call this function
 *            in order to correctly setup the task scheduler.
 */
HEADER_PREFIX void start_scheduler(uint cpu_threads, uint gpu_threads);

/**
 * @brief Stop the task scheduler.
 *
 * @attention This function has to be called by the master and by the slave
 *            thread in order to terminate correctly.
 */
HEADER_PREFIX void stop_scheduler();

/**
 * @brief Macro that simplifies the use of the task scheduler.
 *
 * @param cpu_threads Number of CPU worker threads. This parameter will be passed
 *        directly to @ref start_scheduler.
 * @param gpu_threads Number of CPU threads, that employ the available GPUs.
 *        This parameter will be directly passed to @ref start_scheduler.
 * @param func The function call that shall be handled by the task scheduler.
 *        All necessary parameters of @p func should be passed afterwards.
 */
#define SCHEDULE_OPENCL(cpu_threads, gpu_threads, func, ...) \
    _Pragma("omp parallel num_threads(2)") \
    { \
  start_scheduler(cpu_threads, gpu_threads); \
  if (omp_get_thread_num() == 0) { \
    func(__VA_ARGS__); \
  } \
  stop_scheduler(); \
    }

/**
 * @}
 */

#endif

#endif /* OPENCL_H_ */
