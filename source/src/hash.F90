!--------------------------------------------------------------------------
! Copyright(C) 2008-2016, Kensuke IWAHASHI.  All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!--------------------------------------------------------------------------
!>
!! \file
!! \brief  Module and subroutines to manipulate hashes.
!! \details  This subroutines offer hash manipulation in Fortran90.
!!   To use uses subroutines, you must describe "use hash" in your code.
!! \attention  Currently, hash type is only integer8.
!<
!>
!! \brief Fortran hash module
!<
      module hash
      private
      public maxkey
      public hash_i8
      public hash_new
      public hash_insert
      public hash_lookup
      public hash_delete
      public hash_defined
      public hash_keys
      public hash_destroy
!     integer(4),parameter,private :: primenum=7, primenum2=3
!     integer(4),parameter,private :: primenum=32029, primenum2=37
      integer(4),parameter,private :: primenum=999983, primenum2=37
      integer(4),parameter :: maxkey=128 !< Max length of key
      character(len=256) :: lasterrmsg
      type hash_table_i8  !< A Hash item which has a key and a value.
        character(len=maxkey) :: key  !< Hash key.
        integer(8) :: value  !< Hash value.
        type(hash_table_i8),pointer :: next  !< Array used in conflict.
      end type hash_table_i8
      type hash_i8
        type(hash_table_i8),pointer :: table(:)  !< Hash array table.
      end type hash_i8
      interface hash_new  !< Create new hash.
        module procedure hash_new_i8
      end interface hash_new
      interface hash_insert  !< Add key and value in hash.
        module procedure hash_insert_i8
      end interface hash_insert
      interface hash_delete  !< Remove key and value in hash.
        module procedure hash_delete_i8
      end interface hash_delete
      interface hash_lookup  !< Search value which has specified key.
        module procedure hash_lookup_i8
      end interface hash_lookup
      interface hash_defined  !< Check whether hash key is defined.
        module procedure hash_defined_i8
      end interface hash_defined
      interface hash_keys  !< Return all hash keys.
        module procedure hash_keys_i8
      end interface hash_keys
      interface hash_destroy  !< Destroy hash.
        module procedure hash_destroy_i8
      end interface hash_destroy
      contains
!----------------------------------------------------------------------
!>
!! \brief  Function to calculate hash value from key.
!! \param[in] string  Strin of hash key.
!! \param[out] lerr  True if error, otherwise False.
!! \return  Hash value
!<
      integer(4) function hash_value(string, lerr)
      character(len=*), intent(in) :: string
      logical, intent(out) :: lerr
      integer(8) :: h
      integer(4) :: i, n, h4

      lerr = .false.
      n = len_trim(string)
      if (n > primenum) then
        write(lasterrmsg, '(a, a)')'Too long hash key.  ', string
        lerr = .true.
        return
      endif
      h = 0
      do i=1, n
        h = h * primenum2 + ichar(string(i:i))
      enddo
      h4 = h
      hash_value = mod(abs(h4), primenum)+1

      end function hash_value
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to initialize hash.
!! \param[in,out] a  Hash variable.
!<
      subroutine hash_new_i8(a)
      implicit none
      type(hash_i8) :: a
      integer(4) :: i
      allocate(a%table(primenum))
      do i=1, primenum
        nullify(a%table(i)%next)
      enddo
      end subroutine hash_new_i8
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to insert a key and a value into hash.
!! \param[in] a  Hash variable.
!! \param[in] key  Hash key.
!! \param[in] value  Value corresponding to hash key.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine hash_insert_i8(a, key, value, lerr)
      implicit none
      type(hash_i8) :: a
      character(len=*),intent(in) :: key
      integer(8),intent(in) :: value
      logical, intent(out) :: lerr
      type(hash_table_i8),pointer :: p, new, prev
      integer(4) :: hval

      lerr = .false.
      hval = hash_value(key, lerr)
      prev => a%table(hval)
      p => a%table(hval)%next
      do while (associated(p))
        if (trim(p%key) .eq. trim(key))  exit
        prev => p
        p => p%next
      enddo
      if (associated(p)) then
        p%value = value
      else
        allocate(new)
        new%key = key
        new%value = value
        nullify(new%next)
        prev%next => new
      endif

      end subroutine hash_insert_i8
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to return hash value corresponding to hash key.
!! \param[in] a  Hash variable.
!! \param[in] key  Hash key.
!! \param[out] value  Value corresponding to hash key.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine hash_lookup_i8(a, key, value, lerr)
      implicit none
      type(hash_i8) :: a
      character(len=*),intent(in) :: key
      integer(8),intent(out) :: value
      logical, intent(out) :: lerr
      integer(4) :: hval
      type(hash_table_i8),pointer :: p

      lerr = .false.
      hval = hash_value(key, lerr)
      p => a%table(hval)%next
      do while (associated(p))
        if (trim(p%key) .eq. trim(key))  exit
        p => p%next
      enddo
      if (associated(p)) then
        value = p%value
      else
        lasterrmsg = 'There is no key in hash.'
        value = 0
        lerr = .true.
      endif

      end subroutine hash_lookup_i8
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to remove a key and a value from hash.
!! \param[in,out] a  Hash variable.
!! \param[in] key  Hash key.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine hash_delete_i8(a, key, lerr)
      implicit none
      type(hash_i8) :: a
      character(len=*),intent(in) :: key
      logical, intent(out) :: lerr
      integer(4) :: hval
      type(hash_table_i8),pointer :: p, prev

      lerr = .false.
      hval = hash_value(key, lerr)
      prev => a%table(hval)
      p => a%table(hval)%next
      do while (associated(p))
        if (trim(p%key) .eq. trim(key))  exit
        prev => p
        p => p%next
      enddo
      if (associated(p)) then
        prev%next => p%next
        deallocate(p)
      else
        lasterrmsg = 'There is no key in hash.'
        lerr = .true.
      endif

      end subroutine hash_delete_i8
!----------------------------------------------------------------------
!>
!! \brief  Function to check whether hash key is defined.
!! \param[in] a  Hash variable.
!! \param[in] key  Hash key.
!! \param[out] lerr  True if error, otherwise False.
!<
      logical function hash_defined_i8(a, key, lerr)
      type(hash_i8) :: a
      character(len=*),intent(in) :: key
      logical, intent(out) :: lerr
      integer(4) :: hval
      type(hash_table_i8),pointer :: p
      lerr = .false.
      hval = hash_value(key, lerr)
      p => a%table(hval)%next
      do while (associated(p))
        if (trim(p%key) .eq. trim(key))  exit
        p => p%next
      enddo
      hash_defined_i8 = associated(p)
      end function hash_defined_i8
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to return array of hash keys.
!! \param[in] a  Hash variable.
!! \param[out] keys  Arrays of hash keys.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine hash_keys_i8(a, keys, lerr)
      use array
      implicit none
      type(hash_i8) :: a
      type(array_string), intent(out) :: keys
      type(hash_table_i8),pointer :: p
      logical, intent(out) :: lerr
      integer(4) :: i
      lerr = .false.
      do i=1, primenum
        p => a%table(i)%next
        do while (associated(p))
          call array_push(keys, p%key)
          p => p%next
        enddo
      enddo
      end subroutine hash_keys_i8
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to destroy hash variable.
!! \param[in,out] a  Hash variable.
!<
      subroutine hash_destroy_i8(a)
      implicit none
      type(hash_i8) :: a
      integer(4) :: i
      type(hash_table_i8),pointer :: p, next
      do i=1, primenum
        p => a%table(i)%next
        do while (associated(p))
          next => p%next
          deallocate(p)
          p => next
        enddo
      enddo
      deallocate(a%table)
      end subroutine hash_destroy_i8
!----------------------------------------------------------------------
      end module hash
!======================================================================
!      program m
!      use hash
!      implicit none
!      type(hash_i8) :: a
!      integer(8) :: i1=100, i2
!      call hash_new_i8(a)
!      call hash_insert_i8(a, 'ab1', i1)
!      call hash_insert_i8(a, 'ab8', i1)
!      call hash_lookup_i8(a, 'ab8', i2)
!      write(6,*) i2
!      call hash_destroy_i8(a)
!      end program
