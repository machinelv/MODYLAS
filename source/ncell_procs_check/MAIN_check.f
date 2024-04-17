c======================================================================
      module input_parameters
c======================================================================
      integer(4) :: nprocs
      integer(4) :: nlevel
      integer(4) :: npz,npy,npx
      integer(4) :: ncellx,ncelly,ncellz
      integer(4) :: nxdiv,nydiv,nzdiv
      integer(4) :: lxdiv,lydiv,lzdiv
      logical :: mpi_manual_division_flg=.true.
      character(6) :: dividemethod
      end module
c======================================================================
      program main
c======================================================================
      use input_parameters

      call read_condition()
      call init_cmm_domain_div()

      stop
      end
c----------------------------------------------------------------------
      subroutine read_condition()
c----------------------------------------------------------------------
      use input_parameters

      open(1,file='./CONDITION')
      read(1,*) nlevel
      read(1,*) ncellx,ncelly,ncellz
      read(1,*) nprocs
      read(1,'(a6)') dividemethod
      if(dividemethod(1:4)=='auto')then
        mpi_manual_division_flg=.false.
      else
        read(1,*) nxdiv,nydiv,nzdiv
      endif
      close(1)

      if(nprocs/=nxdiv*nydiv*nzdiv)then
      write(*,*) 'ERROR: nprocs .ne. npx*npy*npz'
      stop
      endif
     
      return
      end
c----------------------------------------------------------------------
      subroutine init_cmm_domain_div()
c----------------------------------------------------------------------
      use input_parameters
      implicit none
      integer(4) :: maxdiv
      integer(4) :: i,j,k
      integer(4) :: mindiv,mxdiv,mydiv,mzdiv

      maxdiv = nprocs
      IF(mpi_manual_division_flg)THEN
        continue
      ELSE
        nxdiv=1
        nydiv=1
        nzdiv=1
      ENDIF
!ya
!     if(nprocs.ne.1) then
!ya
        IF(mpi_manual_division_flg)THEN
          maxdiv=1
          continue
        ELSE
!ya->
          DO WHILE (mod(maxdiv,3).eq.0)
          !!== 3-power ==!!
          nxdiv = nxdiv*3
          maxdiv = maxdiv/3
          if(mod(maxdiv,3).eq.0) then
            nydiv = nydiv*3
            maxdiv = maxdiv/3
          end if
          if(mod(maxdiv,3).eq.0) then
            nzdiv = nzdiv*3
            maxdiv = maxdiv/3
          end if
          END DO 

          !!== 2-power ==!!
          DO WHILE (mod(maxdiv,2).eq.0)
          mindiv=min(nxdiv,nydiv,nzdiv)
          if(    mindiv==nxdiv)then
            nxdiv = nxdiv*2
            maxdiv = maxdiv/2
          elseif(mindiv==nydiv)then
            nydiv = nydiv*2
            maxdiv = maxdiv/2
          else
            nzdiv = nzdiv*2
            maxdiv = maxdiv/2
          endif
          END DO
!ya<-
        ENDIF
!ya
!     end if

!     if(((nxdiv*nydiv*nzdiv).ne.nprocs) .or.
!    &   (maxdiv.ne.1) ) then
!       write(*,*) ' -error init_cmm_domain_div 1-'
!       write(*,*) 'The Number of MPI-procs is incorrect. '
!       write(*,*) 'nxdiv,nydiv,nzdiv,nprocs= ',nxdiv,nydiv,nzdiv,nprocs
!       stop
!     end if

!     maxdiv = max(maxdiv,nxdiv)
!     maxdiv = max(maxdiv,nydiv)
!     maxdiv = max(maxdiv,nzdiv)

!     if(maxdiv.gt.min(ncellx,ncelly,ncellz)) then
!       write(*,*) ' -error init_cmm_domain_div 2-'
!       write(*,*) ' maxdiv.gt.ncell'
!       stop
!     endif

!     mxdiv = mod(ncellx,nxdiv)
!     mydiv = mod(ncelly,nydiv)
!     mzdiv = mod(ncellz,nzdiv)

!     if((mxdiv.ne.0).or.(mydiv.ne.0).or.(mzdiv.ne.0)) then
!       write(*,*) ' -error init_cmm_domain_div 3-'
!       write(*,*) ' mod(ncellxyz,mpidiv).ne.0'
!       stop
!     endif

      lxdiv = ncellx/nxdiv
      lydiv = ncelly/nydiv
      lzdiv = ncellz/nzdiv
!
      IF(mpi_manual_division_flg)THEN
        write(*,*) 'MPI manual division'
      ELSE
        write(*,*) 'MPI auto division'
      ENDIF

      write(*,*) 'nprocs              =', nprocs
      write(*,*) 'ncellx,ncelly,ncellz=', ncellx,ncelly,ncellz
      write(*,*) 'npx   ,npy   ,npz   =', nxdiv,nydiv,nzdiv
      write(*,*) 'lxdiv ,lydiv ,lzdiv =', lxdiv,lydiv,lzdiv

      call  check_supercell_versus_nprocs() !! by Ichikawa

      write(*,*) '====================================='
      write(*,*) 'Check successfully ended.            '
      write(*,*) 'Condition you determined above is OK.'
      write(*,*) '====================================='

      return
      end

c----------------------------------------------------------------------
      subroutine check_supercell_versus_nprocs()
c----------------------------------------------------------------------
!###
!### check subcell and supercell structure ###!
!### by Ichikawa @ 2012/12/26
!###
      use input_parameters
      implicit none
      integer(4) :: i,j,k
      integer(4) :: mx,my,mz,nx,ny,nz,nl,il,iu,jl,ju,kl,ku
      integer(4) :: nscellx,nscelly,nscellz
      integer(4) :: nscxdiv,nscydiv,nsczdiv
      integer(4) :: mcell_size_x,mcell_size_y,mcell_size_z

      npx = nxdiv
      npy = nydiv
      npz = nzdiv
! factorize sub-cell number.
      i = 0; mx = -1
      do while(i == 0)
        mx = mx + 1; i = mod(ncellx/(2**mx),2)
      end do
      i = 0; nx = -1
      do while(i == 0)
        nx = nx + 1; i = mod(ncellx/(2**mx*3**nx),3)
      end do
      j = 0; my = -1
      do while(j == 0)
        my = my + 1; j = mod(ncelly/(2**my),2)
      end do
      j = 0; ny = -1
      do while(j == 0)
        ny = ny + 1; j = mod(ncelly/(2**my*3**ny),3)
      end do
      k = 0; mz = -1
      do while(k == 0)
        mz = mz + 1; k = mod(ncellz/(2**mz),2)
      end do
      k = 0; nz = -1
      do while(k == 0)
        nz = nz + 1; k = mod(ncellz/(2**mz*3**nz),3)
      end do
      if(2**mx*3**nx /= ncellx .or.
     &   2**my*3**ny /= ncelly .or. 2**mz*3**nz /= ncellz) then
        write(6,*) "++ Error. Number of sub-cells is not ",
     &       "the product of power of 2 and power of 3 number. ",
     &       "ncellx: ",ncellx," ncelly: ",ncelly," ncellz: ",ncellz
        stop
      endif

      DO nl = 0, nlevel    ! ---  level -loop.  ---

! cell size in an unit of sub-cell.
      il = nl
      if(il >  mx) il = mx
      iu = nl
      if(iu <= mx) iu = mx
      if(iu > mx + nx) iu = mx + nx
      mcell_size_x = 2**il * 3**(iu - mx)
      jl = nl
      if(jl >  my) jl = my
      ju = nl
      if(ju <= my) ju = my
      if(ju > my + ny) ju = my + ny
      mcell_size_y = 2**jl * 3**(ju - my)
      kl = nl
      if(kl >  mz) kl = mz
      ku = nl
      if(ku <= mz) ku = mz
      if(ku > mz + nz) ku = mz + nz
      mcell_size_z = 2**kl * 3**(ku - mz)

      nscellx = (ncellx - 1) / mcell_size_x + 1
      nscelly = (ncelly - 1) / mcell_size_y + 1
      nscellz = (ncellz - 1) / mcell_size_z + 1
      nscxdiv = (nscellx - 1) / npx + 1
      nscydiv = (nscelly - 1) / npy + 1
      nsczdiv = (nscellz - 1) / npz + 1

      if(nscellx > npx) then
        if(mod(nscellx,npx) /= 0) then
          write(6,*) "++ Error. Wrong process number npx; ",
     &     npx," : super-cells nscellx; ",nscellx," cannot be ",
     &     "divided evenly by processes."
          stop
        endif
      endif
      if(nscellx < npx) then
        if(mod(npx,nscellx) /= 0) then
          write(6,*) "++ Error. Wrong process number npx; ",
     &     npx," : process number cannot be divided evenly by cells ",
     &     "nscellx; ",nscellx,"."
          stop
        endif
      endif
      if(nscelly > npy) then
        if(mod(nscelly,npy) /= 0) then
          write(6,*) "++ Error. Wrong process number npy; ",
     &     npy," : super-cells nscelly; ",nscelly," cannot be ",
     &     "divided evenly by processes."
          stop
        endif
      endif
      if(nscelly < npy) then
        if(mod(npy,nscelly) /= 0) then
          write(6,*) "++ Error. Wrong process number npy; ",
     &     npy," : process number cannot be divided evenly by cells ",
     &     "nscelly; ",nscelly,"."
          stop
        endif
      endif
      if(nscellz > npz) then
        if(mod(nscellz,npz) /= 0) then
          write(6,*) "++ Error. Wrong process number npz; ",
     &     npz," : super-cells nscellz; ",nscellz," cannot be ",
     &     "divided evenly by processes."
          stop
        endif
      endif
      if(nscellz < npz) then
        if(mod(npz,nscellz) /= 0) then
          write(6,*) "++ Error. Wrong process number npz; ",
     &     npz," : process number cannot be divided evenly by cells ",
     &     "nscellz; ",nscellz,"."
          stop
        endif
      endif

      END DO

      return
      end
