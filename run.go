package mpi

//#cgo LDFLAGS: -lmpi
//#include <mpi.h>
import "C"
import "unsafe"

func Init() {
	C.MPI_Init(nil, nil)
}

func Comm_size(comm Comm) int {
	var  size C.int
	C.MPI_Comm_size(comm, &size)
	return int(size)
}

func Comm_rank(comm Comm) int {
	var rank C.int
	C.MPI_Comm_rank(comm, &rank)
	return int(rank)
}


func Allreduce(sendbuf *float64,recvbuf *float64, count int,
	datatype Datatype, op Op, comm Comm){
	C.MPI_Allreduce(unsafe.Pointer(sendbuf), unsafe.Pointer(recvbuf), C.int(count),
		datatype, op, comm)
}
