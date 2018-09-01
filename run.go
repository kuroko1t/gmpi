package mpi

//#cgo LDFLAGS: -lmpi
//#include <mpi.h>
import "C"
import "unsafe"



func Abort(comm COMM, errorcode INT) INT {
	C.MPI_Abort(comm, errorcode)
}
func Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
	int target_rank, MPI_Aint target_disp, int target_count,
	MPI_Datatype target_datatype, MPI_Op op, MPI_Win win) {
}
func Add_error_class(int *errorclass);
func Add_error_code(int errorclass, int *errorcode);
func Add_error_string(int errorcode, const char *string);
func Address(void *location, MPI_Aint *address)
                               __mpi_interface_deprecated__("MPI_Address is superseded by MPI_Get_address in MPI-2.0");
func Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, int recvcount,
                                 MPI_Datatype recvtype, MPI_Comm comm);
func Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, int recvcount,
                                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
func Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, const int recvcounts[],
                                  const int displs[], MPI_Datatype recvtype, MPI_Comm comm);
func Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                  void *recvbuf, const int recvcounts[],
                                  const int displs[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
func Alloc_mem(MPI_Aint size, MPI_Info info,
                                 void *baseptr);
func Allreduce(const void *sendbuf, void *recvbuf, int count,
                                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
func Iallreduce(const void *sendbuf, void *recvbuf, int count,
                                 MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
func Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, MPI_Comm comm);
func Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
func Alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                                 MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                                 const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm);
func Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                                 MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                                 const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
func Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
                                 void *recvbuf, const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[],
                                 MPI_Comm comm);
func Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[], const MPI_Datatype sendtypes[],
                                 void *recvbuf, const int recvcounts[], const int rdispls[], const MPI_Datatype recvtypes[],
                                 MPI_Comm comm, MPI_Request *request);
func Attr_delete(MPI_Comm comm, int keyval)
                                   __mpi_interface_deprecated__("MPI_Attr_delete is superseded by MPI_Comm_delete_attr in MPI-2.0");
func Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag)
                                __mpi_interface_deprecated__("MPI_Attr_get is superseded by MPI_Comm_get_attr in MPI-2.0");
func Attr_put(MPI_Comm comm, int keyval, void *attribute_val)
                                __mpi_interface_deprecated__("MPI_Attr_put is superseded by MPI_Comm_set_attr in MPI-2.0");
func Barrier(MPI_Comm comm);
func Ibarrier(MPI_Comm comm, MPI_Request *request);
func Bcast(void *buffer, int count, MPI_Datatype datatype,
                             int root, MPI_Comm comm);
func Bsend(const void *buf, int count, MPI_Datatype datatype,
                             int dest, int tag, MPI_Comm comm);
func Ibcast(void *buffer, int count, MPI_Datatype datatype,
				                              int root, MPI_Comm comm,
											  MPI_Request *request);
func Bsend_init(const void *buf, int count, MPI_Datatype datatype,
                                  int dest, int tag, MPI_Comm comm, MPI_Request *request);
func Buffer_attach(void *buffer, int size);
func Buffer_detach(void *buffer, int *size);
func Cancel(MPI_Request *request);
func Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]);
func Cart_create(MPI_Comm old_comm, int ndims, const int dims[],
                                   const int periods[], int reorder, MPI_Comm *comm_cart);
func Cart_get(MPI_Comm comm, int maxdims, int dims[],
                                int periods[], int coords[]);
func Cart_map(MPI_Comm comm, int ndims, const int dims[],
                                const int periods[], int *newrank);
func Cart_rank(MPI_Comm comm, const int coords[], int *rank);
func Cart_shift(MPI_Comm comm, int direction, int disp,
                                  int *rank_source, int *rank_dest);
func Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *new_comm);
func Cartdim_get(MPI_Comm comm, int *ndims);
func Close_port(const char *port_name);
func Comm_accept(const char *port_name, MPI_Info info, int root,
                                   MPI_Comm comm, MPI_Comm *newcomm);
func  MPI_Ffunc Comm_c2f(MPI_Comm comm);
func Comm_call_errhandler(MPI_Comm comm, int errorcode);
func Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
func Comm_connect(const char *port_name, MPI_Info info, int root,
                                    MPI_Comm comm, MPI_Comm *newcomm);
func Comm_create_errhandler(MPI_Comm_errhandler_function *function,
                                              MPI_Errhandler *errhandler);
func Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                                          MPI_Comm_delete_attr_function *comm_delete_attr_fn,
                                          int *comm_keyval, void *extra_state);
func Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm);
func Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
func Comm_delete_attr(MPI_Comm comm, int comm_keyval);
func Comm_disconnect(MPI_Comm *comm);
func Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
func Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request);
func Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm);
func  MPI_Comm MPI_Comm_f2c(MPI_Fint comm);
func Comm_free_keyval(int *comm_keyval);
func Comm_free(MPI_Comm *comm);
func Comm_get_attr(MPI_Comm comm, int comm_keyval,
                                     void *attribute_val, int *flag);
func Dist_graph_create(MPI_Comm comm_old, int n, const int nodes[],
                                         const int degrees[], const int targets[],
                                         const int weights[], MPI_Info info,
                                         int reorder, MPI_Comm * newcomm);
func Dist_graph_create_adjacent(MPI_Comm comm_old,
                                                  int indegree, const int sources[],
                                                  const int sourceweights[],
                                                  int outdegree,
                                                  const int destinations[],
                                                  const int destweights[],
                                                  MPI_Info info, int reorder,
                                                  MPI_Comm *comm_dist_graph);
func func Dist_graph_neighbors(MPI_Comm comm, int maxindegree,
                                           int sources[], int sourceweights[],
                                           int maxoutdegree,
                                           int destinations[],
                                           int destweights[]);
func Dist_graph_neighbors_count(MPI_Comm comm,
                                                  int *inneighbors,
                                                  int *outneighbors,
                                                  int *weighted);
func Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *erhandler);
func Comm_get_info(MPI_Comm comm, MPI_Info *info_used);
func Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen);
func Comm_get_parent(MPI_Comm *parent);
func Comm_group(MPI_Comm comm, MPI_Group *group);
func Comm_join(int fd, MPI_Comm *intercomm);
func Comm_rank(MPI_Comm comm, int *rank);
func Comm_remote_group(MPI_Comm comm, MPI_Group *group);
func Comm_remote_size(MPI_Comm comm, int *size);
func Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val);
func Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
func Comm_set_info(MPI_Comm comm, MPI_Info info);
func Comm_set_name(MPI_Comm comm, const char *comm_name);
func Comm_size(MPI_Comm comm, int *size);
func Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info,
                                  int root, MPI_Comm comm, MPI_Comm *intercomm,
                                  int array_of_errcodes[]);
func Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                                           const int array_of_maxprocs[], const MPI_Info array_of_info[],
                                           int root, MPI_Comm comm, MPI_Comm *intercomm,
                                           int array_of_errcodes[]);
func Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
func Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
func Comm_test_inter(MPI_Comm comm, int *flag);
func Compare_and_swap(const void *origin_addr, const void *compare_addr,
                                        void *result_addr, MPI_Datatype datatype, int target_rank,
                                        MPI_Aint target_disp, MPI_Win win);
func Dims_create(int nnodes, int ndims, int dims[]);
func  MPI_Ffunc Errhandler_c2f(MPI_Errhandler errhandler);
func Errhandler_create(MPI_Handler_function *function,
                                         MPI_Errhandler *errhandler)
                                         __mpi_interface_deprecated__("MPI_Errhandler_create is superseded by MPI_Comm_create_errhandler in MPI-2.0");
func  MPI_Errhandler MPI_Errhandler_f2c(MPI_Fint errhandler);
func Errhandler_free(MPI_Errhandler *errhandler);
func Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
                                      __mpi_interface_deprecated__("MPI_Errhandler_get is superseded by MPI_Comm_get_errhandler in MPI-2.0");
func Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
                                      __mpi_interface_deprecated__("MPI_Errhandler_set is superseded by MPI_Comm_set_errhandler in MPI-2.0");
func Error_class(int errorcode, int *errorclass);
func Error_string(int errorcode, char *string, int *resultlen);
func Exscan(const void *sendbuf, void *recvbuf, int count,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
func Fetch_and_op(const void *origin_addr, void *result_addr, MPI_Datatype datatype,
                                    int target_rank, MPI_Aint target_disp, MPI_Op op, MPI_Win win);
func Iexscan(const void *sendbuf, void *recvbuf, int count,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
func  MPI_Ffunc File_c2f(MPI_File file);
func  MPI_File MPI_File_f2c(MPI_Fint file);
func File_call_errhandler(MPI_File fh, int errorcode);
func File_create_errhandler(MPI_File_errhandler_function *function,
                                              MPI_Errhandler *errhandler);
func File_set_errhandler( MPI_File file, MPI_Errhandler errhandler);
func File_get_errhandler( MPI_File file, MPI_Errhandler *errhandler);
func File_open(MPI_Comm comm, const char *filename, int amode,
                                 MPI_Info info, MPI_File *fh);
func File_close(MPI_File *fh);
func File_delete(const char *filename, MPI_Info info);
func File_set_size(MPI_File fh, MPI_Offset size);
func File_preallocate(MPI_File fh, MPI_Offset size);
func File_get_size(MPI_File fh, MPI_Offset *size);
func File_get_group(MPI_File fh, MPI_Group *group);
func File_get_amode(MPI_File fh, int *amode);
func File_set_info(MPI_File fh, MPI_Info info);
func File_get_info(MPI_File fh, MPI_Info *info_used);
func File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,
                                     MPI_Datatype filetype, const char *datarep, MPI_Info info);
func File_get_view(MPI_File fh, MPI_Offset *disp,
                                     MPI_Datatype *etype,
                                     MPI_Datatype *filetype, char *datarep);
func File_read_at(MPI_File fh, MPI_Offset offset, void *buf,
                                    int count, MPI_Datatype datatype, MPI_Status *status);
func File_read_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                                        int count, MPI_Datatype datatype, MPI_Status *status);
func File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,
                                     int count, MPI_Datatype datatype, MPI_Status *status);
func File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
                                         int count, MPI_Datatype datatype, MPI_Status *status);
func File_iread_at(MPI_File fh, MPI_Offset offset, void *buf,
                                     int count, MPI_Datatype datatype, MPI_Request *request);
func File_iwrite_at(MPI_File fh, MPI_Offset offset, const void *buf,
                                      int count, MPI_Datatype datatype, MPI_Request *request);
func File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                                     int count, MPI_Datatype datatype, MPI_Request *request);
func File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf,
                                      int count, MPI_Datatype datatype, MPI_Request *request);
func File_read(MPI_File fh, void *buf, int count,
                                 MPI_Datatype datatype, MPI_Status *status);
func File_read_all(MPI_File fh, void *buf, int count,
                                     MPI_Datatype datatype, MPI_Status *status);
func File_write(MPI_File fh, const void *buf, int count,
                                  MPI_Datatype datatype, MPI_Status *status);
func File_write_all(MPI_File fh, const void *buf, int count,
                                      MPI_Datatype datatype, MPI_Status *status);
func File_iread(MPI_File fh, void *buf, int count,
                                  MPI_Datatype datatype, MPI_Request *request);
func File_iwrite(MPI_File fh, const void *buf, int count,
                                   MPI_Datatype datatype, MPI_Request *request);
func File_iread_all(MPI_File fh, void *buf, int count,
                                  MPI_Datatype datatype, MPI_Request *request);
func File_iwrite_all(MPI_File fh, const void *buf, int count,
                                   MPI_Datatype datatype, MPI_Request *request);
func File_seek(MPI_File fh, MPI_Offset offset, int whence);
func File_get_position(MPI_File fh, MPI_Offset *offset);
func File_get_byte_offset(MPI_File fh, MPI_Offset offset,
                                            MPI_Offset *disp);
func File_read_shared(MPI_File fh, void *buf, int count,
                                        MPI_Datatype datatype, MPI_Status *status);
func File_write_shared(MPI_File fh, const void *buf, int count,
					 MPI_Datatype datatype, MPI_Status *status);
func File_iread_shared(MPI_File fh, void *buf, int count,
                                         MPI_Datatype datatype, MPI_Request *request);
func File_iwrite_shared(MPI_File fh, const void *buf, int count,
                                          MPI_Datatype datatype, MPI_Request *request);
func File_read_ordered(MPI_File fh, void *buf, int count,
                                         MPI_Datatype datatype, MPI_Status *status);
func File_write_ordered(MPI_File fh, const void *buf, int count,
                                          MPI_Datatype datatype, MPI_Status *status);
func File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
func File_get_position_shared(MPI_File fh, MPI_Offset *offset);
func File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
                                              int count, MPI_Datatype datatype);
func File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
func File_write_at_all_begin(MPI_File fh, MPI_Offset offset, const void *buf,
                                               int count, MPI_Datatype datatype);
func File_write_at_all_end(MPI_File fh, const void *buf, MPI_Status *status);
func File_read_all_begin(MPI_File fh, void *buf, int count,
                                           MPI_Datatype datatype);
func File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
func File_write_all_begin(MPI_File fh, const void *buf, int count,
                                            MPI_Datatype datatype);
func File_write_all_end(MPI_File fh, const void *buf, MPI_Status *status);
func File_read_ordered_begin(MPI_File fh, void *buf, int count,
                                               MPI_Datatype datatype);
func File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
func File_write_ordered_begin(MPI_File fh, const void *buf, int count,
                                                MPI_Datatype datatype);
func File_write_ordered_end(MPI_File fh, const void *buf, MPI_Status *status);
func File_get_type_extent(MPI_File fh, MPI_Datatype datatype,
                                            MPI_Aint *extent);
func File_set_atomicity(MPI_File fh, int flag);
func File_get_atomicity(MPI_File fh, int *flag);
func File_sync(MPI_File fh);
func Finalize(void);
func Finalized(int *flag);
func Free_mem(void *base);
func Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, int recvcount, MPI_Datatype recvtype,
                              int root, MPI_Comm comm);
func Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, int recvcount, MPI_Datatype recvtype,
                              int root, MPI_Comm comm, MPI_Request *request);
func Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const int recvcounts[], const int displs[],
                               MPI_Datatype recvtype, int root, MPI_Comm comm);
func Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, const int recvcounts[], const int displs[],
                               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
func Get_address(const void *location, MPI_Aint *address);
func Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
func Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count);
func Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count);
func Get(void *origin_addr, int origin_count,
                           MPI_Datatype origin_datatype, int target_rank,
                           MPI_Aint target_disp, int target_count,
                           MPI_Datatype target_datatype, MPI_Win win);
func Get_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                                      void *result_addr, int result_count, MPI_Datatype result_datatype,
                                      int target_rank, MPI_Aint target_disp, int target_count,
                                      MPI_Datatype target_datatype, MPI_Op op, MPI_Win win);
func Get_library_version(char *version, int *resultlen);
func Get_processor_name(char *name, int *resultlen);
func Get_version(int *version, int *subversion);
func Graph_create(MPI_Comm comm_old, int nnodes, const int index[],
                                    const int edges[], int reorder, MPI_Comm *comm_graph);
func Graph_get(MPI_Comm comm, int maxindex, int maxedges,
                                 int index[], int edges[]);
func Graph_map(MPI_Comm comm, int nnodes, const int index[], const int edges[],
                                 int *newrank);
func Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
func Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors,
                                       int neighbors[]);
func Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
func Grequest_complete(MPI_Request request);
func Grequest_start(MPI_Grequest_query_function *query_fn,
                                      MPI_Grequest_free_function *free_fn,
                                      MPI_Grequest_cancel_function *cancel_fn,
                                      void *extra_state, MPI_Request *request);
func  MPI_Ffunc Group_c2f(MPI_Group group);
func Group_compare(MPI_Group group1, MPI_Group group2, int *result);
func Group_difference(MPI_Group group1, MPI_Group group2,
                                        MPI_Group *newgroup);
func Group_excl(MPI_Group group, int n, const int ranks[],
                                  MPI_Group *newgroup);
func  MPI_Group MPI_Group_f2c(MPI_Fint group);
func Group_free(MPI_Group *group);
func Group_incl(MPI_Group group, int n, const int ranks[],
                                  MPI_Group *newgroup);
func Group_intersection(MPI_Group group1, MPI_Group group2,
                                          MPI_Group *newgroup);
func Group_range_excl(MPI_Group group, int n, int ranges[][3],
                                        MPI_Group *newgroup);
func Group_range_incl(MPI_Group group, int n, int ranges[][3],
                                        MPI_Group *newgroup);
func Group_rank(MPI_Group group, int *rank);
func Group_size(MPI_Group group, int *size);
func Group_translate_ranks(MPI_Group group1, int n, const int ranks1[],
                                             MPI_Group group2, int ranks2[]);
func Group_union(MPI_Group group1, MPI_Group group2,
                                   MPI_Group *newgroup);
func Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest,
                              int tag, MPI_Comm comm, MPI_Request *request);
func Improbe(int source, int tag, MPI_Comm comm,
                               int *flag, MPI_Message *message,
                               MPI_Status *status);
func Imrecv(void *buf, int count, MPI_Datatype type,
                              MPI_Message *message, MPI_Request *request);
func  MPI_Ffunc Info_c2f(MPI_Info info);
func Info_create(MPI_Info *info);
func Info_delete(MPI_Info info, const char *key);
func Info_dup(MPI_Info info, MPI_Info *newinfo);
func  MPI_Info MPI_Info_f2c(MPI_Fint info);
func Info_free(MPI_Info *info);
func Info_get(MPI_Info info, const char *key, int valuelen,
                                char *value, int *flag);
func Info_get_nkeys(MPI_Info info, int *nkeys);
func Info_get_nthkey(MPI_Info info, int n, char *key);
func Info_get_valuelen(MPI_Info info, const char *key, int *valuelen,
                                         int *flag);
func Info_set(MPI_Info info, const char *key, const char *value);
func Init(int *argc, char ***argv);
func Initialized(int *flag);
func Init_thread(int *argc, char ***argv, int required,
                                   int *provided);
func Intercomm_create(MPI_Comm local_comm, int local_leader,
                                        MPI_Comm bridge_comm, int remote_leader,
                                        int tag, MPI_Comm *newintercomm);
func Intercomm_merge(MPI_Comm intercomm, int high,
                                       MPI_Comm *newintercomm);
func Iprobe(int source, int tag, MPI_Comm comm, int *flag,
                              MPI_Status *status);
func Irecv(void *buf, int count, MPI_Datatype datatype, int source,
                             int tag, MPI_Comm comm, MPI_Request *request);
func Irsend(const void *buf, int count, MPI_Datatype datatype, int dest,
                              int tag, MPI_Comm comm, MPI_Request *request);
func Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
                             int tag, MPI_Comm comm, MPI_Request *request);
func Issend(const void *buf, int count, MPI_Datatype datatype, int dest,
                              int tag, MPI_Comm comm, MPI_Request *request);
func Is_thread_main(int *flag);
func Keyval_create(MPI_Copy_function *copy_fn,
                                     MPI_Delete_function *delete_fn,
                                     int *keyval, void *extra_state)
                                     __mpi_interface_deprecated__("MPI_Keyval_create is superseded by MPI_Comm_create_keyval in MPI-2.0");
func Keyval_free(int *keyval)
                                   __mpi_interface_deprecated__("MPI_Keyval_free is superseded by MPI_Comm_free_keyval in MPI-2.0");
func Lookup_name(const char *service_name, MPI_Info info, char *port_name);
func  MPI_Ffunc Message_c2f(MPI_Message message);
func  MPI_Message MPI_Message_f2c(MPI_Fint message);
func Mprobe(int source, int tag, MPI_Comm comm,
                               MPI_Message *message,
                               MPI_Status *status);
func Mrecv(void *buf, int count, MPI_Datatype type,
                             MPI_Message *message, MPI_Status *status);
func Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                          MPI_Comm comm);
func Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                           void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                           MPI_Comm comm, MPI_Request *request);
func Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                           void *recvbuf, const int recvcounts[], const int displs[],
                                           MPI_Datatype recvtype, MPI_Comm comm);
func Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                            void *recvbuf, const int recvcounts[], const int displs[],
                                            MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
func Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                         void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                         MPI_Comm comm);
func Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                          MPI_Comm comm, MPI_Request *request);
func Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],  MPI_Datatype sendtype,
                                          void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                                          MPI_Comm comm);
func Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                                           void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                                           MPI_Comm comm, MPI_Request *request);
func Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                          void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                                          MPI_Comm comm);
func Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                                           void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                                           MPI_Comm comm, MPI_Request *request);
func  MPI_Ffunc Op_c2f(MPI_Op op);
func Op_commutative(MPI_Op op, int *commute);
func Op_create(MPI_User_function *function, int commute, MPI_Op *op);
func Open_port(MPI_Info info, char *port_name);
func  MPI_Op MPI_Op_f2c(MPI_Fint op);
func Op_free(MPI_Op *op);
func Pack_external(const char datarep[], const void *inbuf, int incount,
                                     MPI_Datatype datatype, void *outbuf,
                                     MPI_Aint outsize, MPI_Aint *position);
func Pack_external_size(const char datarep[], int incount,
                                          MPI_Datatype datatype, MPI_Aint *size);
func Pack(const void *inbuf, int incount, MPI_Datatype datatype,
                            void *outbuf, int outsize, int *position, MPI_Comm comm);
func Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm,
                                 int *size);
func Pcontrol(const int level, ...);
func Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
func Publish_name(const char *service_name, MPI_Info info,
                                    const char *port_name);
func Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                           int target_rank, MPI_Aint target_disp, int target_count,
                           MPI_Datatype target_datatype, MPI_Win win);
func Query_thread(int *provided);
func Raccumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                                   int target_rank, MPI_Aint target_disp, int target_count,
                                   MPI_Datatype target_datatype, MPI_Op op, MPI_Win win, MPI_Request *request);
func Recv_init(void *buf, int count, MPI_Datatype datatype, int source,
                                 int tag, MPI_Comm comm, MPI_Request *request);
func Recv(void *buf, int count, MPI_Datatype datatype, int source,
                            int tag, MPI_Comm comm, MPI_Status *status);
func Reduce(const void *sendbuf, void *recvbuf, int count,
                              MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
func Ireduce(const void *sendbuf, void *recvbuf, int count,
                              MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request *request);
func Reduce_local(const void *inbuf, void *inoutbuf, int count,
                                    MPI_Datatype datatype, MPI_Op op);
func Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
func Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
func Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
func Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
func Register_datarep(const char *datarep,
                                        MPI_Datarep_conversion_function *read_conversion_fn,
                                        MPI_Datarep_conversion_function *write_conversion_fn,
                                        MPI_Datarep_extent_function *dtype_file_extent_fn,
                                        void *extra_state);
func  MPI_Ffunc Request_c2f(MPI_Request request);
func  MPI_Request MPI_Request_f2c(MPI_Fint request);
func Request_free(MPI_Request *request);
func Request_get_status(MPI_Request request, int *flag,
                                          MPI_Status *status);
func Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                            int target_rank, MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype,
                            MPI_Win win, MPI_Request *request);
func Rget_accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                                       void *result_addr, int result_count, MPI_Datatype result_datatype,
                                       int target_rank, MPI_Aint target_disp, int target_count,
                                       MPI_Datatype target_datatype, MPI_Op op,
                                       MPI_Win win, MPI_Request *request);
func Rput(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                            int target_rank, MPI_Aint target_disp, int target_cout,
                            MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request);
func Rsend(const void *ibuf, int count, MPI_Datatype datatype, int dest,
                             int tag, MPI_Comm comm);
func Rsend_init(const void *buf, int count, MPI_Datatype datatype,
                                  int dest, int tag, MPI_Comm comm,
                                  MPI_Request *request);
func Scan(const void *sendbuf, void *recvbuf, int count,
                            MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
func Iscan(const void *sendbuf, void *recvbuf, int count,
                            MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
func Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype,
                               int root, MPI_Comm comm);
func Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype,
                               int root, MPI_Comm comm, MPI_Request *request);
func Scatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                                MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, int root, MPI_Comm comm);
func Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                                MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request);
func Send_init(const void *buf, int count, MPI_Datatype datatype,
                                 int dest, int tag, MPI_Comm comm,
                                 MPI_Request *request);
func Send(const void *buf, int count, MPI_Datatype datatype, int dest,
                            int tag, MPI_Comm comm);
func Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                int dest, int sendtag, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, int source, int recvtag,
                                MPI_Comm comm,  MPI_Status *status);
func Sendrecv_replace(void * buf, int count, MPI_Datatype datatype,
                                        int dest, int sendtag, int source, int recvtag,
                                        MPI_Comm comm, MPI_Status *status);
func Ssend_init(const void *buf, int count, MPI_Datatype datatype,
                                  int dest, int tag, MPI_Comm comm,
                                  MPI_Request *request);
func Ssend(const void *buf, int count, MPI_Datatype datatype, int dest,
                             int tag, MPI_Comm comm);
func Start(MPI_Request *request);
func Startall(int count, MPI_Request array_of_requests[]);
func Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status);
func Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status);
func Status_set_cancelled(MPI_Status *status, int flag);
func Status_set_elements(MPI_Status *status, MPI_Datatype datatype,
                                           int count);
func Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype,
                                             MPI_Count count);
func Testall(int count, MPI_Request array_of_requests[], int *flag,
                               MPI_Status array_of_statuses[]);
func Testany(int count, MPI_Request array_of_requests[], int *index,
                               int *flag, MPI_Status *status);
func Test(MPI_Request *request, int *flag, MPI_Status *status);
func Test_cancelled(const MPI_Status *status, int *flag);
func Testsome(int incount, MPI_Request array_of_requests[],
                                int *outcount, int array_of_indices[],
                                MPI_Status array_of_statuses[]);
func Topo_test(MPI_Comm comm, int *status);
func  MPI_Ffunc Type_c2f(MPI_Datatype datatype);
func Type_commit(MPI_Datatype *type);
func Type_contiguous(int count, MPI_Datatype oldtype,
                                       MPI_Datatype *newtype);
func Type_create_darray(int size, int rank, int ndims,
                                          const int gsize_array[], const int distrib_array[],
                                          const int darg_array[], const int psize_array[],
                                          int order, MPI_Datatype oldtype,
                                          MPI_Datatype *newtype);
func Type_create_f90_complex(int p, int r, MPI_Datatype *newtype);
func Type_create_f90_integer(int r, MPI_Datatype *newtype);
func Type_create_f90_real(int p, int r, MPI_Datatype *newtype);
func Type_create_hindexed_block(int count, int blocklength,
                                                  const MPI_Aint array_of_displacements[],
                                                  MPI_Datatype oldtype,
                                                  MPI_Datatype *newtype);
func Type_create_hindexed(int count, const int array_of_blocklengths[],
                                            const MPI_Aint array_of_displacements[],
                                            MPI_Datatype oldtype,
                                            MPI_Datatype *newtype);
func Type_create_hvector(int count, int blocklength, MPI_Aint stride,
                                           MPI_Datatype oldtype,
                                           MPI_Datatype *newtype);
func Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                                          MPI_Type_delete_attr_function *type_delete_attr_fn,
                                          int *type_keyval, void *extra_state);
func Type_create_indexed_block(int count, int blocklength,
                                                 const int array_of_displacements[],
                                                 MPI_Datatype oldtype,
                                                 MPI_Datatype *newtype);
func Type_create_struct(int count, const int array_of_block_lengths[],
                                          const MPI_Aint array_of_displacements[],
                                          const MPI_Datatype array_of_types[],
                                          MPI_Datatype *newtype);
func Type_create_subarray(int ndims, const int size_array[], const int subsize_array[],
                                            const int start_array[], int order,
                                            MPI_Datatype oldtype, MPI_Datatype *newtype);
func Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb,
                                           MPI_Aint extent, MPI_Datatype *newtype);
func Type_delete_attr(MPI_Datatype type, int type_keyval);
func Type_dup(MPI_Datatype type, MPI_Datatype *newtype);
func Type_extent(MPI_Datatype type, MPI_Aint *extent)
                                   __mpi_interface_deprecated__("MPI_Type_extent is superseded by MPI_Type_get_extent in MPI-2.0");
func Type_free(MPI_Datatype *type);
func Type_free_keyval(int *type_keyval);
func  MPI_Datatype MPI_Type_f2c(MPI_Fint datatype);
func Type_get_attr(MPI_Datatype type, int type_keyval,
                                     void *attribute_val, int *flag);
func Type_get_contents(MPI_Datatype mtype, int max_integers,
                                         int max_addresses, int max_datatypes,
                                         int array_of_integers[],
                                         MPI_Aint array_of_addresses[],
                                         MPI_Datatype array_of_datatypes[]);
func Type_get_envelope(MPI_Datatype type, int *num_integers,
                                         int *num_addresses, int *num_datatypes,
                                         int *combiner);
func Type_get_extent(MPI_Datatype type, MPI_Aint *lb,
                                       MPI_Aint *extent);
func Type_get_extent_x(MPI_Datatype type, MPI_Count *lb,
                                         MPI_Count *extent);
func Type_get_name(MPI_Datatype type, char *type_name,
                                     int *resultlen);
func Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb,
                                            MPI_Aint *true_extent);
func Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *true_lb,
                                              MPI_Count *true_extent);
func Type_hindexed(int count, int array_of_blocklengths[],
                                     MPI_Aint array_of_displacements[],
                                     MPI_Datatype oldtype, MPI_Datatype *newtype)
                                     __mpi_interface_deprecated__("MPI_Type_hindexed is superseded by MPI_Type_create_hindexed in MPI-2.0");
func Type_hvector(int count, int blocklength, MPI_Aint stride,
                                    MPI_Datatype oldtype, MPI_Datatype *newtype)
                                    __mpi_interface_deprecated__("MPI_Type_hvector is superseded by MPI_Type_create_hvector in MPI-2.0");
func Type_indexed(int count, const int array_of_blocklengths[],
                                    const int array_of_displacements[],
                                    MPI_Datatype oldtype, MPI_Datatype *newtype);
func Type_lb(MPI_Datatype type, MPI_Aint *lb)
                               __mpi_interface_deprecated__("MPI_Type_lb is deprecated, use MPI_Type_get_extent in MPI-2.0");
func Type_match_size(int typeclass, int size, MPI_Datatype *type);
func Type_set_attr(MPI_Datatype type, int type_keyval,
                                     void *attr_val);
func Type_set_name(MPI_Datatype type, const char *type_name);
func Type_size(MPI_Datatype type, int *size);
func Type_size_x(MPI_Datatype type, MPI_Count *size);
func Type_struct(int count, int array_of_blocklengths[],
                                   MPI_Aint array_of_displacements[],
                                   MPI_Datatype array_of_types[],
                                   MPI_Datatype *newtype)
                                   __mpi_interface_deprecated__("MPI_Type_struct is superseded by MPI_Type_create_struct in MPI-2.0");
func Type_ub(MPI_Datatype mtype, MPI_Aint *ub)
                               __mpi_interface_deprecated__("MPI_Type_ub is deprecated, use MPI_Type_get_extent in MPI-2.0");
func Type_vector(int count, int blocklength, int stride,
                                   MPI_Datatype oldtype, MPI_Datatype *newtype);
func Unpack(const void *inbuf, int insize, int *position,
                              void *outbuf, int outcount, MPI_Datatype datatype,
                              MPI_Comm comm);
func Unpublish_name(const char *service_name, MPI_Info info, const char *port_name);
func Unpack_external (const char datarep[], const void *inbuf, MPI_Aint insize,
                                        MPI_Aint *position, void *outbuf, int outcount,
                                        MPI_Datatype datatype);
func Waitall(int count, MPI_Request array_of_requests[],
                               MPI_Status *array_of_statuses);
func Waitany(int count, MPI_Request array_of_requests[],
                               int *index, MPI_Status *status);
func Wait(MPI_Request *request, MPI_Status *status);
func Waitsome(int incount, MPI_Request array_of_requests[],
                                int *outcount, int array_of_indices[],
                                MPI_Status array_of_statuses[]);
func Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                                    MPI_Comm comm, void *baseptr, MPI_Win *win);
func Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info,
                                           MPI_Comm comm, void *baseptr, MPI_Win *win);
func Win_attach(MPI_Win win, void *base, MPI_Aint size);
func  MPI_Ffunc Win_c2f(MPI_Win win);
func Win_call_errhandler(MPI_Win win, int errorcode);
func Win_complete(MPI_Win win);
func Win_create(void *base, MPI_Aint size, int disp_unit,
                                  MPI_Info info, MPI_Comm comm, MPI_Win *win);
func Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win);
func Win_create_errhandler(MPI_Win_errhandler_function *function,
                                             MPI_Errhandler *errhandler);
func Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                                         MPI_Win_delete_attr_function *win_delete_attr_fn,
                                         int *win_keyval, void *extra_state);
func Win_delete_attr(MPI_Win win, int win_keyval);
func Win_detach(MPI_Win win, const void *base);
func  MPI_Win MPI_Win_f2c(MPI_Fint win);
func Win_fence(int assert, MPI_Win win);
func Win_flush(int rank, MPI_Win win);
func Win_flush_all(MPI_Win win);
func Win_flush_local(int rank, MPI_Win win);
func Win_flush_local_all(MPI_Win win);
func Win_free(MPI_Win *win);
func Win_free_keyval(int *win_keyval);
func Win_get_attr(MPI_Win win, int win_keyval,
                                    void *attribute_val, int *flag);
func Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler);
func Win_get_group(MPI_Win win, MPI_Group *group);
func Win_get_info(MPI_Win win, MPI_Info *info_used);
func Win_get_name(MPI_Win win, char *win_name, int *resultlen);
func Win_lock(int lock_type, int rank, int assert, MPI_Win win);
func Win_lock_all(int assert, MPI_Win win);
func Win_post(MPI_Group group, int assert, MPI_Win win);
func Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val);
func Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler);
func Win_set_info(MPI_Win win, MPI_Info info);
func Win_set_name(MPI_Win win, const char *win_name);
func Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
func Win_start(MPI_Group group, int assert, MPI_Win win);
func Win_sync(MPI_Win win);
func Win_test(MPI_Win win, int *flag);
func Win_unlock(int rank, MPI_Win win);
func Win_unlock_all(MPI_Win win);
func Win_wait(MPI_Win win);
func  double MPI_Wtick(void);
func  double MPI_Wtime(void);
