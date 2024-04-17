/*----------------------------------------------------------------------
*MODYLAS ver. 1.1.0 
*
*Copyright (c) 2014-2019 Nagoya University
*              2020-2023 The University of Tokyo
*
*Released under the MIT license.
*see https://opensource.org/licenses/MIT
*-----------------------------------------------------------------------
*MODYLAS Developers:
*Yoshimichi Andoh, Kazushi Fujimoto, Tatsuya Sakashita, Noriyuki Yoshii, 
*Zhiye Tang, Jiachao Zhang, Yuta Asano, Ryo Urano, Tetsuro Nagai, 
*Atsushi Yamada, Hidekazu Kojima, Kensuke Iwahashi, Fumiyasu Mizutani, 
*Shin-ichi Ichikawa, and Susumu Okazaki.
*----------------------------------------------------------------------*/
#include <mpi-ext.h>
#include <stdio.h>
#include <stdlib.h>
#define RDMA_INTERVAL          (8192)
#define HEAP_SIZE              (512*1024*1024) // 512 MB
#define RDMA_TAG               (0)
#define RDMA_FLAG              (FJMPI_RDMA_LOCAL_NIC0 | FJMPI_RDMA_REMOTE_NIC0 | FJMPI_RDMA_PATH0)
#define ALIGN                  (32)
#define MEMID_TABLE_SIZE       (509)
#define HEAP_OFFSET_TABLE_SIZE (7)
#define DUMMY_VAL              (0)
static int _put_num = 0;
static uint64_t _post_laddr, *_post_raddr;
static int _post_buf, _post_memid = 510, *_wait_num;
static void *_heap_buf = NULL;
static int _heap_memid = 509;
static uint64_t *_heap_raddr;
static uint64_t _heap_laddr;
static uint64_t _heap_offset = 0;
static uint64_t _heap_offset_table[HEAP_OFFSET_TABLE_SIZE];
static int _heap_offset_num = 0;
static struct FJMPI_Rdma_cq _cq;
typedef struct element{
  void     *key_addr;
  uint64_t laddr;
  uint64_t *raddr;
} element_t;
static element_t _table[MEMID_TABLE_SIZE];
static int _memid = 0;
static int _size, _rank;
static MPI_Comm _comm;

void rdma_get_remote_addr_(const int memid, uint64_t *raddr){
  // Wait until rdma_reg_mem_..() is completed on other ranks 
  // before calling this function.
  MPI_Barrier(_comm);

  // Although FJMPI_Rdma_get_remote_addr() is called using a while() statement
  // in the document of FJRDMA, the algorithm makes the network heavy.
  // Thus, we use the following light algorithm.
  for(int ncount=0,i=1; i<_size+1; ncount++,i++){
    int target = (_rank + _size - i) % _size;
    raddr[target] = FJMPI_Rdma_get_remote_addr(target, memid);
    if(raddr[target] == FJMPI_RDMA_ERROR){
      printf("Error : Getting remote address of RDMA.\n");
      MPI_Finalize();
      exit(1);
    }
    if(ncount > RDMA_INTERVAL){
      MPI_Barrier(_comm);
      ncount = 0;
    }
  }
}

void rdma_put_(const int rank, const uint64_t raddr, const uint64_t laddr, const size_t length){
  FJMPI_Rdma_put(rank, RDMA_TAG, raddr, laddr, length, RDMA_FLAG);
  _put_num++;
}

void rdma_put_post_(const int rank, const uint64_t raddr, const uint64_t laddr, const size_t length){
  FJMPI_Rdma_put(rank, RDMA_TAG, raddr, laddr, length, RDMA_FLAG | FJMPI_RDMA_REMOTE_NOTICE);
  _put_num++;
}

static void rdma_sync_memory(){
  while(1){
    if(_put_num == 0) break;
    int tmp = FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC0, &_cq);
    if(tmp == FJMPI_RDMA_NOTICE) _put_num--;
    else if(tmp == FJMPI_RDMA_HALFWAY_NOTICE) _wait_num[_cq.pid]--;
  }
}

void rdma_sync_all_(){
  rdma_sync_memory();
  MPI_Barrier(_comm);
}

void rdma_post_(const int rank){
  size_t length = sizeof(int);
  FJMPI_Rdma_put(rank, RDMA_TAG, _post_raddr[rank], _post_laddr, length, RDMA_FLAG | FJMPI_RDMA_REMOTE_NOTICE);
  _put_num++;
}

void rdma_wait_(const int rank){
  _wait_num[rank]++;
  while(1){
    if(_wait_num[rank] <= 0) break;
    int tmp = FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC0, &_cq);
    if(tmp == FJMPI_RDMA_NOTICE) _put_num--;
    else if(tmp == FJMPI_RDMA_HALFWAY_NOTICE) _wait_num[_cq.pid]--;
  }
  rdma_sync_memory();
}

void rdma_init_(const MPI_Fint fcomm){
  FJMPI_Rdma_init();

  // Prepare POST/WAIT
  _comm =  MPI_Comm_f2c(fcomm);
  MPI_Comm_size(_comm, &_size);
  MPI_Comm_rank(_comm, &_rank);
  _post_raddr = malloc(sizeof(uint64_t) * _size);
  _post_laddr = FJMPI_Rdma_reg_mem(_post_memid, &_post_buf, sizeof(int));
  rdma_get_remote_addr_(_post_memid, _post_raddr);

  _heap_buf   = malloc(HEAP_SIZE);
  _heap_raddr = malloc(sizeof(uint64_t) * _size);
  _heap_laddr = FJMPI_Rdma_reg_mem(_heap_memid, _heap_buf, HEAP_SIZE);
  rdma_get_remote_addr_(_heap_memid, _heap_raddr);

  _wait_num = malloc(sizeof(int) * _size);
  for(int i=0;i<_size;i++)
    _wait_num[i] = 0;
}

void rdma_finalize_(){
  FJMPI_Rdma_finalize();
}

static void rdma_register_addr( void *addr, const size_t length)
{
  if(_memid >= MEMID_TABLE_SIZE){
    if(_rank == 0)
      fprintf(stderr, "_memid is over\n", addr);
    MPI_Abort(_comm, 1);
  }
  _table[_memid].key_addr = addr;
  _table[_memid].laddr    = FJMPI_Rdma_reg_mem(_memid, addr, length);
  _table[_memid].raddr    = malloc(sizeof(uint64_t) * _size);
  rdma_get_remote_addr_(_memid, _table[_memid].raddr);
  _memid++;
}

void rdma_register_addr_i1_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}

void rdma_register_addr_i3_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}

void rdma_register_addr_r2_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}

void rdma_register_addr_c1_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}

void rdma_register_addr_c4_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}

#ifdef SEP_RI
void rdma_register_addr_s5_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}

void rdma_register_addr_r5_(void *addr, const size_t length)
{
  rdma_register_addr(addr, length);
}
#endif

static uint64_t* rdma_get_raddr(const void *addr)
{
  for(int i=0;i<_memid;i++)
    if(_table[i].key_addr == addr)
      return _table[i].raddr;

  if(_rank == 0)
    fprintf(stderr, "key (%p) is not registered\n", addr);

  MPI_Abort(_comm, 1);
  return DUMMY_VAL;
}

uint64_t* rdma_get_raddr_i1_(const void *addr)
{
  return rdma_get_raddr(addr);
}

uint64_t* rdma_get_raddr_i3_(const void *addr)
{
  return rdma_get_raddr(addr);
}

uint64_t* rdma_get_raddr_r2_(const void *addr)
{
  return rdma_get_raddr(addr);
}

uint64_t* rdma_get_raddr_c1_(const void *addr)
{
  return rdma_get_raddr(addr);
}

uint64_t* rdma_get_raddr_c4_(const void *addr)
{
  return rdma_get_raddr(addr);
}

#ifdef SEP_RI
uint64_t* rdma_get_raddr_s5_(const void *addr)
{
  return rdma_get_raddr(addr);
}

uint64_t* rdma_get_raddr_r5_(const void *addr)
{
  return rdma_get_raddr(addr);
}
#endif

static uint64_t rdma_get_laddr(const void *addr)
{
  for(int i=0;i<_memid;i++)
    if(_table[i].key_addr == addr)
      return _table[i].laddr;

  if(_rank == 0)
    fprintf(stderr, "key (%p) is not registered\n", addr);

  MPI_Abort(_comm, 1);
  return DUMMY_VAL;
}

uint64_t rdma_get_laddr_i1_(const void *addr)
{
  return rdma_get_laddr(addr);
}

uint64_t rdma_get_laddr_i3_(const void *addr)
{
  return rdma_get_laddr(addr);
}

uint64_t rdma_get_laddr_r2_(const void *addr)
{
  return rdma_get_laddr(addr);
}

uint64_t rdma_get_laddr_c1_(const void *addr)
{
 return rdma_get_laddr(addr);
}

uint64_t rdma_get_laddr_c4_(const void *addr)
{
  return rdma_get_laddr(addr);
}

#ifdef SEP_RI
uint64_t rdma_get_laddr_s5_(const void *addr)
{
  return rdma_get_laddr(addr);
}

uint64_t rdma_get_laddr_r5_(const void *addr)
{
  return rdma_get_laddr(addr);
}
#endif

void* rdma_allocate_coarray_(const size_t size, uint64_t *ret_offset)
{
  if(_heap_offset_num >= HEAP_OFFSET_TABLE_SIZE){
    if(_rank == 0)
      fprintf(stderr, "HEAP_OFFSET_TABLE_SIZE is too small.\n");
    MPI_Abort(_comm, 1);
  }

  *ret_offset             = _heap_offset;
  uint64_t current_offset = _heap_offset;
  uint64_t new_offset     = (((_heap_offset+size)+ALIGN-1)/ALIGN)*ALIGN;
  _heap_offset_table[_heap_offset_num++] = new_offset - current_offset;
  _heap_offset = new_offset;

  if(_heap_offset > HEAP_SIZE){
    if(_rank == 0)
      fprintf(stderr, "HEAP_SIZE is too small %d %d\n", 
	      (int)_heap_offset, HEAP_SIZE);
    MPI_Abort(_comm, 1);
  }

  return (char*)_heap_buf + current_offset;
}

void rdma_deallocate_(const int *num)
{
  if(_heap_offset_num < *num){
    if(_rank == 0)
      fprintf(stderr, "Intger value num is too large %d %d\n", 
	      _heap_offset_num, *num);
    MPI_Abort(_comm, 1);
  }

  for(int i=0;i<*num;i++)
    _heap_offset -= _heap_offset_table[--_heap_offset_num];
}

uint64_t rdma_get_heap_laddr_()
{
  return _heap_laddr;
}

uint64_t* rdma_get_heap_raddr_()
{
  return _heap_raddr;
}
