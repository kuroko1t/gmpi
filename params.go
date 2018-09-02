package mpi

//#cgo LDFLAGS: -lmpi
//#include <mpi.h>
import "C"

type COMM C.MPI_Comm
type Aint           C.MPI_Aint
type Offset         C.MPI_Offset
type Count          C.MPI_Count
type Comm           C.MPI_Comm
type Datatype       C.MPI_Datatype
type Errhandler     C.MPI_Errhandler
type File           C.MPI_File
type Group          C.MPI_Group
type Info           C.MPI_Info
type Op             C.MPI_Op
type Request        C.MPI_Request
type Message        C.MPI_Message
type Status         C.MPI_Status
type Win            C.MPI_Win
type T_enum         C.MPI_T_enum
type T_cvar_handle  C.MPI_T_cvar_handle
type T_pvar_handle  C.MPI_T_pvar_handle
type T_pvar_session C.MPI_T_pvar_session

//comm
var COMM_WORLD =  C.MPI_COMM_WORLD
var SUM = C.MPI_SUM
var Float64 = C.MPI_DOUBLE
