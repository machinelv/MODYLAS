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
!! \brief Module and subroutines to manipulate arrays
!<
!>
!! \brief Module to manipulate arrays
!<
      module array
      private
      public array_string
      public array_new
      public array_number
      public array_pop
      public array_push
      public array_join
      public array_nth
      public array_destroy
      type array_string  !< Number and pointers of array_item_string
        integer(4) :: n=0  !< Number of array
        type(array_item_string), pointer:: top  !< Top of array_item_string
        type(array_item_string), pointer:: last !< Last of array_item_string
      end type array_string
      type array_item_string !< Doubly linked list of string.
        character(len=128) :: string !< String value
        type(array_item_string), pointer :: prev !< Previous pointer
        type(array_item_string), pointer :: next !< Next pointer
      end type array_item_string
      interface array_new
        module procedure array_new_string
      end interface array_new
      interface array_number
        module procedure array_number_string
      end interface array_number
      interface array_pop
        module procedure array_pop_string
      end interface array_pop
      interface array_push
        module procedure array_push_string
      end interface array_push
      interface array_join
        module procedure array_join_string
      end interface array_join
      interface array_nth
        module procedure array_nth_string
      end interface array_nth
      interface array_destroy
        module procedure array_destroy_string
      end interface array_destroy
      contains
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to create new string array.
!! \param[out] a  Newly initialized string array.
!<
      subroutine array_new_string(a)
      implicit none
      type(array_string) :: a
      a%n = 0
      nullify(a%top)
      nullify(a%last)
      end subroutine
!----------------------------------------------------------------------
!>
!! \brief  Function to get number of string array length.
!! \param[in] a  String array.
!<
      integer(4) function array_number_string(a)
      implicit none
      type(array_string) :: a
      array_number_string = a%n
      end function array_number_string
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to add string to string array.
!! \param[in,out] a  String array.
!! \param[in] string String to be added.
!<
      subroutine array_push_string(a, string)
      implicit none
      type(array_string) :: a
      character(len=*) :: string
      type(array_item_string),pointer :: new, oldlast
      allocate(new)
      if (a%n == 0) then
        a%top => new
        a%last => new
        nullify(new%prev)
      else
        oldlast => a%last
        oldlast%next => new
        a%last => new
        new%prev => oldlast
      endif
      nullify(new%next)
      new%string = string
      a%n = a%n + 1
!     write(6,*) 'array push: ', a%n, trim(string)
      return
      end subroutine array_push_string
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to get nth string from string arrays.
!! \param[in] a  String array.
!! \param[in] nth  Number of order.
!! \param[out] string  String to be stored.
!<
      subroutine array_nth_string(a, nth, string)
      implicit none
      type(array_string) :: a
      integer(4) :: nth
      character(len=*) :: string
      type(array_item_string),pointer :: p
      integer(4) :: i
      if (nth > a%n .or. nth <= 0) then
        write(6,*) 'Illegal subscription.'
        string = ''
      else
        p => a%top
        do i = 1, nth-1
          p => p%next
        enddo
        string = p%string
      endif
      end subroutine array_nth_string
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to take string from string array.
!! \param[in] a  String array.
!! \param[out] string  String to be stored.
!<
      subroutine array_pop_string(a, string)
      implicit none
      type(array_string) :: a
      character(len=*) :: string
      type(array_item_string),pointer :: oldlast, oldlast2

      if (a%n > 1) then
        string = a%last%string
        oldlast2 => a%last%prev
        oldlast => a%last
        a%last => oldlast2
        nullify(oldlast2%next)
        deallocate(oldlast)
        a%n = a%n - 1
      else if (a%n == 1) then
        string = a%last%string
        nullify(a%top)
        nullify(a%last)
        a%n = 0
      else
        write(6,*) 'array is empty.'
        string = ''
      endif
      
!     write(6,*) 'array pop: ', a%n, trim(string)
      return
      end subroutine array_pop_string
!----------------------------------------------------------------------
!>
!! \brief Subroutine to make string from string array with separator.
!! \param[in] sep  String of separator.
!! \param[in] a  String array.
!! \param[out] s  String to be stored.
!<
      subroutine array_join_string(sep, a, s)
      implicit none
      character(len=*),intent(in) :: sep
      type(array_string) :: a
      character(len=*),intent(out) :: s
      character(len=128) :: item
      integer(4) :: i, n
      type(array_item_string),pointer :: p

      n = array_number_string(a)
      s = ''
      if (a%n > 0) then
        s = a%top%string
        p => a%top%next
        do i = 2, a%n
          s = trim(s) // trim(sep) // p%string
          p => p%next
        enddo
      endif
!     write(6,*) trim(s)

      return
      end subroutine array_join_string
!----------------------------------------------------------------------
!>
!! \brief  Subroutine to destroy string array.
!! \param[in,out] a  String array.
!<
      subroutine array_destroy_string(a)
      implicit none
      type(array_string) :: a
      type(array_item_string),pointer :: oldlast, oldlast2
      if (a%n > 0) then
        oldlast2 => a%last%prev
        oldlast => a%last
        a%last => oldlast2
        deallocate(oldlast)
        a%n = a%n - 1
      endif
      return
      end subroutine array_destroy_string
!----------------------------------------------------------------------
      endmodule array
!======================================================================
!     program array_test
!     use array
!     implicit none
!     type(array_string) :: a
!     character(len=128) :: string
!     call array_new_string(a)
!     call array_push_string(a, 'abc')
!     call array_push_string(a, 'def')
!     call array_push_string(a, 'ghi')
!     write(6,*) array_number_string(a)
!     call array_nth_string(a, 2, string)
!     write(6,*) trim(string)
!     call array_pop_string(a, string)
!     write(6,*) trim(string)
!     call array_pop_string(a, string)
!     write(6,*) trim(string)
!     call array_pop_string(a, string)
!     write(6,*) trim(string)
!     call array_pop_string(a, string)
!     write(6,*) trim(string)
!     call array_destroy_string(a)
!     end program
