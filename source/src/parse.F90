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
!! \brief  Module and subroutines to read input file
!     in Fortran90 flexibly.
!<
!>
!! \brief  Module to read input file in Fortran90 flexibly.
!<
      module parse
      use array
      use hash
      implicit none
      private
      public parse_ignore_key
      public parse_open
#ifdef GROEXT
      public parse_open_string
#endif
      public parse_get
      public parse_specify_key
      public parse_close
      public parse_show_error_message
      public parse_set_debug
      public maxkey
      integer(4),parameter:: ionum=99     !< Temporal unit number.
      character, pointer:: allstr(:)  !< All letters in input.
      character(len=maxkey) :: lastkey !< Key with asterisk and number.
      character(len=256) :: lasterrmsg     !< Last error message.
      integer(8) :: lallstr            !< Location of last of the file.
      type(hash_i8) :: vsta, vend, ignored_key
      logical :: increment
      logical:: use_ignored_key=.false., mode_ignored
      type(array_string) :: rawkeywords
      integer(4):: debug=0, i
      logical:: normal_letter(0:127) = (/(.false.,i=0,32), &
     &                                   (.true., i=33,59), &
     &                                   (.false.,i=60,62), &
     &                                   (.true., i=63,126), &
     &                                    .false./)
      !>
      !! \brief  Interface to get value corresponding to key.
      !<
      interface parse_get
        module procedure parse_get_string
        module procedure parse_get_integer4
        module procedure parse_get_integer8
        module procedure parse_get_real4
        module procedure parse_get_real8
        module procedure parse_get_logical
      end interface parse_get
      interface parse_set_debug
        module procedure parse_set_debug_logical
        module procedure parse_set_debug_integer4
      end interface parse_set_debug
      interface int2string
        module procedure int4_to_string
        module procedure int8_to_string
      end interface int2string
      interface string2int
        module procedure string_to_int4
      end interface string2int
      contains
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to read and parse input file.
!! \param[in] filename  Input file.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_open(filename, lerr)
      use mpi_tool
      implicit none
      character(LEN=*), intent(in) :: filename
      logical, intent(out) :: lerr
      integer(4) :: istatus
      character(len=128) :: s
      character(len=9):: form_type
#ifdef MPI
      include 'mpif.h'
      integer(4) :: ierr, myrank
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
#endif  /* MPI */
      lastkey = '/'
      lallstr = -1
      lerr = .false.
      increment=.false.
      mode_ignored=.false.
      call array_new(rawkeywords)
      call array_push(rawkeywords, '')
#ifdef MPI
      if (myrank == 0) then
#endif  /* MPI */
        lasterrmsg = ''
#ifdef __GFORTRAN__
        form_type = 'formatted'
#else /* __GFORTRAN__ */
        form_type = 'binary'
#endif /* __GFORTRAN__ */
        open(ionum, file=filename, status='old', form=form_type, iostat=istatus)
        if (istatus /= 0) then
          write(0, '(a,a,a)') 'ERROR: Cannot open for reading  ', filename
          call modylas_abort()
        endif
        call parse_open__load
        close(ionum)
        call array_pop(rawkeywords, s)
        if (array_number(rawkeywords) > 0) then
          call array_join('/', rawkeywords, s)
          write(lasterrmsg, '(a)') 'Tag is not closed.  ' // s
          lerr = .true.
          return
        endif
#ifdef MPI
      endif
      call MPI_BCAST(lallstr,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(allstr, lallstr, MPI_CHARACTER, &
     &               0,MPI_COMM_WORLD,ierr)
#endif  /* MPI */
      call parse_analyze(lerr)
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (Internal subroutine) Count number of valid letters in the
!!    specified unit.
!<
      subroutine parse_open__load
      use array
      implicit none
#ifdef __GFORTRAN__
      integer(4), parameter:: size_line=1024
#else /* __GFORTRAN__ */
      integer(4), parameter:: size_line=1024
#endif /* __GFORTRAN__ */
      integer(4) :: istatus, si, so, i, n_output, ic, nc, j
      integer(8):: lsize,lsizenew
      character, pointer :: allstrnew(:)
      character :: c
      character(len=128):: output
      character(len=size_line):: line
      logical :: ignore
      istatus = 0
      si = 0
      lallstr = 0
      ignore = .false.
      lsize = 10000
      allocate(allstr(lsize))
      do while (istatus == 0 .or. istatus == -2)
        line = ''
#ifdef __GFORTRAN__
        read(ionum, '(a)', advance='no', iostat=istatus, &
     &        size=nc) line(2:size_line)
!       Add space character at the begining of the line.
        nc = nc + 1
        line(1:1) = ' '
        do ic = 1, nc
#else /* __GFORTRAN__ */
        read(ionum, iostat=istatus) line
        do ic = 1, size_line
#endif /* __GFORTRAN__ */
          c = line(ic:ic)
          if (ichar(c) == 10) then
            ignore = .false.
          endif
          if (ignore) then
            cycle
          endif
          if (ichar(c) <= 32) then
            c = ' '
          endif
          if (c == '#') then
            ignore = .true.
            cycle
          else if (c == ';' .or. c == ',') then
            c = ' '
          endif
          call parse_open__state_value(c, si, so, output, n_output)
          if (lallstr + n_output > lsize) then
            lsizenew = lsize + lsize
            do while (lallstr + n_output > lsizenew)
              lsizenew = lsizenew + lsizenew
            enddo
            allocate(allstrnew(lsizenew))
!           bug using intel compiler
!            allstrnew(:lsize) = allstr(:lsize)
            do j = 1, lsize
              allstrnew(j) = allstr(j)
            enddo 
            deallocate(allstr)
            allstr => allstrnew
            lsize = lsizenew
          endif
          do i = 1, n_output
            allstr(lallstr+i) = output(i:i)
          enddo
          lallstr = lallstr + n_output
          si = so
        enddo
#ifdef __GFORTRAN__
        ignore = .false.
#endif /* __GFORTRAN__ */
      enddo
      if (debug >= 2) then
        call parse_write_allstr(allstr, lallstr)
      endif
      end subroutine parse_open__load
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (internal subroutine) Return an input state.
!! \param[in] c  Current letter to parse.
!! \param[in] si  Last state before parse.
!! \param[out] so  Result of state after parse.
!! \param[out] output  Decided string.
!! \param[out] n_output Length of decided string.
!! \details
!!   Table of state is following.
!!   \code
!!
!!  Consider 16 patterns in inputs.
!!          <a>  </a> val  k=v
!!   <a>     o     o    o    o
!!   </a>    o     o    o    o
!!   val     o     o    o    o
!!   k=v     o     o    o    o
!!
!!                     so
!!          0 10 11 12 13 20 21 22 23 99
!!      0  ss kv -- -- -- k- -- -- --  o    (Empty)
!!     10  -- ss kv K- -- kV -- -- --  o    a
!!     11  -- -V ss K- -- kV -- -- --  o    a (space)
!!     12  -- -- -- ss -v -- -- -- --  o    a=
!! si  13  -V -- -- -- ss kV -- -- --  o    a=A
!!     20  -- -- -- -- -- -- k- k- --  o    <
!!     21  -K -- -- -- -- -- ss -- --  o    <a
!!     22  -- -- -- -- -- -- -- -- k-  o    </
!!     23  -K -- -- -- -- -- -- -- ss  o    </a
!!
!!   k: Not decided keyword
!!   K: Decided keyword
!!   v: Not decided value
!!   V: Decided value
!!  ss: Recurrence
!!   \endcode
!<
      subroutine parse_open__state_value(c, si, so, output, n_output)
      use array
      implicit none
      character, intent(in)::c
      integer(4), intent(in):: si
      integer(4), intent(out):: so
      character(len=128), intent(out):: output
      integer(4), intent(out):: n_output
      character(len=72),save:: str=''
      character(len=128):: s
      integer(4),save:: length_str=0
      logical:: keep, lerr
!
      keep = .false.
      output = ''
      n_output = 0
      if (si == 0) then
        if (c == ' ') then
          so = 0
          return
        else if (normal_letter(ichar(c))) then
          so = 10
          keep = .true.
        else if (c == '<') then
          so = 20
        else
          so = 99
        endif
      else if (si == 10) then
        if (c == ' ') then
          so = 11
          keep = .true.
        else if (c == '=') then
          so = 12
!         Decide keyword
          call array_push(rawkeywords, str)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          if (debug >= 1) then
            write(6,*) trim(s), mode_ignored
          endif
          if (mode_ignored) then
            output = ''
            n_output = 0
            keep = .true.
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (c == '<') then
          so = 20
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (normal_letter(ichar(c))) then
          so = 10
          keep = .true.
        else
          so = 99
        endif
      else if (si == 11) then
        if (c == ' ') then
          so = 11
          return
        else if (c == '=') then
          so = 12
!         -- Decide keyword --
          call array_push(rawkeywords, str)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          if (debug >= 1) then
            write(6,*) trim(s), mode_ignored
          endif
          if (mode_ignored) then
            output = ''
            n_output = 0
            keep = .true.
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (c == '<') then
          so = 20
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (normal_letter(ichar(c))) then
          so = 10
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
          keep = .true.
        else
          so = 99
        endif
      else if (si == 12) then
        if (c == ' ') then
          so = 12
          return
        else if (normal_letter(ichar(c))) then
          so = 13
          keep = .true.
        else
          so = 99
        endif
      else if (si == 13) then
        if (c == ' ') then
          so = 0
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
          call array_pop(rawkeywords, s)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
        else if (normal_letter(ichar(c))) then
          so = 13
          keep = .true.
        else
          so = 99
        endif
      else if (si == 20) then
        if (c == '/') then
          so = 22
        else
          so = 21
          keep = .true.
        endif
      else if (si == 21) then
        if (c == '>') then
          so = 0
!         -- Decide keyword --
          call array_push(rawkeywords, str)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          if (debug >= 1) then
            write(6,*) trim(s), mode_ignored
          endif
          output = str
          n_output = length_str
          str = ''
          length_str = 0
        else
          so = 21
          keep = .true.
        endif
      else if (si == 22) then
        if (c == '>') then
          so = 99
        else
          so = 23
          keep = .true.
        endif
      else if (si == 23) then
        if (c == '>') then
          so = 0
!         -- Decide keyword --
          call array_pop(rawkeywords, s)
          if (s(:length_str) .ne. str(:length_str)) then
            write(6,*) 'Warning: Keyword differ.'
            write(6,*) ' Start: "' // trim(s) // '".'
            write(6,*) ' End: "' // str(:length_str) // '".'
          endif
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          output = str
          n_output = length_str
          str = ''
          length_str = 0
        else
          so = 23
          keep = .true.
        endif
      endif
      if (keep) then
        length_str = length_str + 1
        str(length_str:length_str) = c
      else
        n_output = n_output + 1
        output(n_output:n_output) = c
      endif
      if (debug >= 2) then
        write(6,*) c, si, so, '"', str(:length_str), '" "', &
     &             output(:n_output), '"'
      endif
      end subroutine parse_open__state_value
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (Intenal subroutine) Parse tags.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_analyze(lerr)
      use array
      use hash
      implicit none
      logical, intent(out) :: lerr
      character(len=128) :: s,keyword,merged,lastkeyword
      character(len=8) :: snum
      logical :: flag_equal=.false.
      type(array_string) :: indexedkeywords
      type(hash_i8) :: keyvarray
      integer(8) :: m, msta, mend, ksta, kend, is, ie, p, peq
      lerr = .true.
      call array_new(indexedkeywords)
      call array_push(indexedkeywords, '')
      call hash_new(keyvarray)
      call hash_new(vsta)
      call hash_new(vend)
      is = 1
      ie = lallstr
      do while (is > 0)
        is = parse_verify(allstr, is, ie, ' ')
        if (is == 0) then
!         case : End with space.
          exit
        endif
        if (allstr(is) == '<') then
          if (is+1 > ie) then
            lasterrmsg = 'End with <.'
            goto 900
          endif
          msta = is
          mend = parse_scan(allstr, is, ie, '>')
          if (mend == 0) then
            lasterrmsg = 'Parenthesis error.'
            goto 900
          endif
          if (allstr(msta+1) == '/') then
!           case : </keyword>
            ksta = msta + 2
            kend = mend - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            if (flag_equal) then
              call array_pop(indexedkeywords, s)
              call array_pop(rawkeywords, s)
              flag_equal = .false.
            endif
            call array_pop(indexedkeywords, s)
            call array_pop(rawkeywords, s)
            if (trim(s) /= trim(keyword)) then
              write(lasterrmsg, '(a,a,a,a)') &
     &            'Parenthesis keyword differ.', &
     &            trim(s), ' ', trim(keyword)
              goto 900
            endif
          else
!           case : <keyword>
            ksta = msta + 1
            kend = mend - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            if (flag_equal) then
              call array_pop(indexedkeywords, s)
              call array_pop(rawkeywords, s)
              flag_equal = .false.
            endif
            call array_join('/', indexedkeywords, merged)
            if (ichar(merged(2:2)) < 32) then
              merged = '/' // keyword
            else
              merged = trim(merged) // '/' // keyword
            endif
            call hash_lookup(keyvarray, merged, m, lerr)
            m = m + 1
            call hash_insert(keyvarray, merged, m, lerr)
            call int2string(m, snum)
            s = trim(keyword) // '*' // snum
            call array_push(indexedkeywords, s)
            call array_push(rawkeywords, keyword)
            lastkeyword = keyword
          endif
          is = mend + 1
        else
          p = parse_scan(allstr, is, ie, '<>= ')
          if (allstr(p) == ' ') then
!           Remove extra spaces before equal
            peq = parse_verify(allstr, p, ie, ' ')
            if (allstr(peq) == '=') then
              p = peq
            endif
          endif
          if (p == 0) then
!           case : End with value string.
            ksta = is
            kend = ie
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            is = 0
          else if (allstr(p) == '>') then
            lasterrmsg = 'Value contains >.'
            goto 900
          else if (allstr(p) == '<' .or. allstr(p) == ' ') then
!           case : value
            ksta = is
            kend = p - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            is = p
            call array_pop(indexedkeywords, s)
            call array_pop(rawkeywords, s)
            call array_join('/', indexedkeywords, merged)
            if (ichar(merged(2:2)) < 32) then
              merged = '/' // lastkeyword
            else
              merged = trim(merged) // '/' // lastkeyword
            endif
            call hash_lookup(keyvarray, merged, m, lerr)
            call hash_insert(keyvarray, merged, m+1, lerr)
            call int2string(m, snum)
            s = trim(lastkeyword) // '*' // snum
            call array_push(indexedkeywords, s)
            call array_push(rawkeywords, lastkeyword)
            call array_join('/', indexedkeywords, s)
            call hash_insert(vsta, s, ksta, lerr)
            call hash_insert(vend, s, kend, lerr)
          else if (allstr(p) == '=') then
!           case : keyword=
            ksta = is
            kend = p - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            is = p + 1
            if (flag_equal) then
              call array_pop(indexedkeywords, s)
              call array_pop(rawkeywords, s)
              flag_equal = .false.
            endif
            call array_join('/', indexedkeywords, merged)
            if (ichar(merged(2:2)) < 32) then
              merged = '/' // keyword
            else
              merged = trim(merged) // '/' // keyword
            endif
            call hash_lookup(keyvarray, merged, m, lerr)
            m = m + 1
            call hash_insert(keyvarray, merged, m, lerr)
            call int2string(m, snum)
            s = trim(keyword) // '*' // snum
            call array_push(indexedkeywords, s)
            call array_push(rawkeywords, keyword)
            lastkeyword = keyword
            flag_equal = .true.
          endif
        endif
      enddo
      call array_pop(indexedkeywords, s)
      if (array_number(indexedkeywords) /= 0) then
        lasterrmsg = 'ERROR: Unexpected EOF.'
        goto 900
      endif
      lerr = .false.
  900 continue
      call array_destroy(indexedkeywords)
      call array_destroy(rawkeywords)
      call hash_destroy(keyvarray)
      return
      end subroutine parse_analyze
      end subroutine parse_open
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to write partial string in debug.
!! \param[in] a  Partial string array to write.
!! \param[in] lsize  Length to be written.
!<
      subroutine parse_write_allstr(a,lsize)
      implicit none
      character, intent(in) :: a(:)
      integer(8), intent(in):: lsize
      write(6,*) a(:lsize)
      end subroutine parse_write_allstr
! ---------------------------------------------------------------------
!>
!! \brief  (Correspond to 'scan' in Fortran90)
!! \param[in] a  String array.
!! \param[in] is0  Location of head of string array.
!! \param[in] ie0  Location of tail of string array.
!! \param[in] substr  A set of characters.
!! \param[in] back  Flag to scan from tail.
!! \return  First location of presence.
!<
      integer(8) function parse_scan(a, is0, ie0, substr, back)
      implicit none
      character, intent(in) :: a(:)
      integer(8), intent(in) :: is0, ie0
      character(LEN=*), intent(in) :: substr
      logical, intent(in), optional :: back
      integer(8) :: i, is, ie
      integer(4) :: j, istep
      character :: c
      parse_scan = 0
      is = is0
      ie = ie0
      istep = 1
      if (present(back)) then
        if (back) then
          is = ie0
          ie = is0
          istep = -1
        endif
      endif
      do i=is, ie, istep
        do j=1,len(substr)
          c = substr(j:j)
          if (a(i) == c) then
            parse_scan = i
            return
          endif
        enddo
      enddo
      end function parse_scan
! ---------------------------------------------------------------------
!>
!! \brief  (Correspond to 'verify' in Fortran90)
!! \param[in] a  String array.
!! \param[in] is0  Location of head of string array.
!! \param[in] ie0  Location of tail of string array.
!! \param[in] substr  A set of characters.
!! \param[in] back  Flag to verify from tail.
!! \return  First location of absence.
!<
      integer(8) function parse_verify(a, is0, ie0, substr, back)
      implicit none
      character, intent(in) :: a(:)
      integer(8), intent(in) :: is0, ie0
      character(LEN=*), intent(in) :: substr
      logical, intent(in), optional :: back
      integer(8) :: i, is, ie
      integer(4) :: j, istep
      character :: c
      is = is0
      ie = ie0
      istep = 1
      parse_verify = 0
      if (present(back)) then
        if (back) then
          is = ie0
          ie = is0
          istep = -1
        endif
      endif
      do i=is, ie, istep
        parse_verify = i
        do j=1,len(substr)
          c = substr(j:j)
          if (a(i) == c) then
            parse_verify = 0
            exit
          endif
        enddo
        if (parse_verify /= 0) then
          exit
        endif
      enddo
      end function parse_verify
! ---------------------------------------------------------------------
!>
!! \brief  (Subroutine to convert character to internal string array.)
!<
      subroutine chars2string(ista, iend, string, lerr)
      implicit none
      integer(8), intent(in) :: ista, iend
      character(len=*), intent(out) :: string
      logical, intent(out) :: lerr
      integer(4) :: i, l, maxlen
      lerr = .false.
      maxlen = len(string)
      do i=1, maxlen
        string(i:i) = ' '
      enddo
      l = iend - ista + 1
      if (l > maxlen) then
        lasterrmsg = 'Too long items.'
        lerr = .true.
        return
      endif
      do i=1, l
        string(i:i) = allstr(ista+i-1)
      enddo
      end subroutine chars2string
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to convert 4-byte integer to characters.
!! \param[in] n  4-byte integer
!! \param[out] string  Characters converted from integer.
      subroutine int4_to_string(n, string)
      use mpi_tool
      implicit none
      integer(4), intent(in):: n
      character(len=*), intent(out):: string
      character:: ct(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
      integer(4):: i, k=0, nr, n1
      if (n < 10) then
        k = 1
      else if (n < 100) then
        k = 2
      else if (n < 1000) then
        k = 3
      else if (n < 10000) then
        k = 4
      else if (n < 100000) then
        k = 5
      else if (n < 1000000) then
        k = 6
      else if (n < 10000000) then
        k = 7
      else if (n < 100000000) then
        k = 8
      else
        write(0,*) 'ERROR: Number is too big. ', n
        call modylas_abort()
      endif
      nr = n
      string = ''
      do i=k,1,-1
        n1 = mod(nr,10)
        string(i:i) = ct(n1)
        nr = (nr - n1) / 10
      enddo
      end subroutine int4_to_string
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to convert 8-byte integer to characters.
!! \param[in] n  8-byte integer
!! \param[out] string  Characters converted from integer.
      subroutine int8_to_string(n, string)
      use mpi_tool
      implicit none
      integer(8), intent(in):: n
      character(len=*), intent(out):: string
      character:: ct(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
      integer(4):: i, k=0, nr, n1
      if (n < 10) then
        k = 1
      else if (n < 100) then
        k = 2
      else if (n < 1000) then
        k = 3
      else if (n < 10000) then
        k = 4
      else if (n < 100000) then
        k = 5
      else if (n < 1000000) then
        k = 6
      else if (n < 10000000) then
        k = 7
      else if (n < 100000000) then
        k = 8
      else
        write(0,*) 'ERROR: Number is too big. ', n
        call modylas_abort()
      endif
      nr = n
      string = ''
      do i=k,1,-1
        n1 = mod(nr,10)
        string(i:i) = ct(n1)
        nr = (nr - n1) / 10
      enddo
      end subroutine int8_to_string
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to convert characters to 4-byte integer.
!! \param[in] string  characters to be converted.
!! \param[out] n  4-byte integer converted from characters.
      subroutine string_to_int4(string, n)
      implicit none
      character(len=*), intent(in):: string
      integer(4), intent(out):: n
      integer(4):: i, ic, ic0
      character:: c
      n = 0
      ic0 = ichar('0')
      do i=1, len(string)
        c = string(i:i)
        if (c < '0' .or. c > '9') then
          exit
        endif
        n = n * 10 + ichar(c) - ic0
      enddo
      end subroutine string_to_int4
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to store value corresponding to key for integer8.
!! \param[in] key  Key
!! \param[out] iret  Integer8 to be stored
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_integer8(key, iret, lerr)
      implicit none
      character(len=*), intent(in) :: key
      integer(8), intent(out) :: iret
      logical, intent(out) :: lerr
      character(len=80) :: buf
      call parse_get_string(key, buf, lerr)
      if (lerr) then
        iret = 0
      else
        read(buf,*) iret
        if (debug >= 1) write(6,*) 'integer8: ', iret
      endif
      end subroutine parse_get_integer8
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to store value corresponding to key for integer4.
!! \param[in] key  Key
!! \param[out] iret  Integer4 to be stored
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_integer4(key, iret, lerr)
      implicit none
      character(len=*), intent(in) :: key
      integer(4), intent(out) :: iret
      logical, intent(out) :: lerr
      character(len=80) :: buf
      call parse_get_string(key, buf, lerr)
      if (lerr) then
        iret = 0
      else
        read(buf,*) iret
        if (debug >= 1) write(6,*) 'integer4: ', iret
      endif
      end subroutine parse_get_integer4
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to store value corresponding to key for real8.
!! \param[in] key  Key
!! \param[out] ret  Real8 to be stored
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_real8(key, ret, lerr)
      implicit none
      character(len=*), intent(in) :: key
      real(8), intent(out) :: ret
      logical, intent(out) :: lerr
      character(len=80) :: buf
      call parse_get_string(key, buf, lerr)
      if (lerr) then
        ret = 0.0d0
      else
        read(buf,*) ret
        if (debug >= 1) write(6,*) 'real8: ', ret
      endif
      end subroutine parse_get_real8
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to store value corresponding to key for real4.
!! \param[in] key  Key
!! \param[out] ret  Real4 to be stored
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_real4(key, ret, lerr)
      implicit none
      character(len=*), intent(in) :: key
      real(4), intent(out) :: ret
      logical, intent(out) :: lerr
      character(len=80) :: buf
      call parse_get_string(key, buf, lerr)
      if (lerr) then
        ret = 0.0
      else
        read(buf,*) ret
        if (debug >= 1) write(6,*) 'real4: ', ret
      endif
      end subroutine parse_get_real4
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to store value corresponding to key for logical.
!! \param[in] key  Key
!! \param[out] lret  Logical value to be stored
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_logical(key, lret, lerr)
      implicit none
      character(len=*),intent(in) :: key
      logical, intent(out) :: lret, lerr
      character(len=80) :: buf
      call parse_get_string(key, buf, lerr)
      if (lerr) then
        lret = .false.
      else
        read(buf,*) lret
        if (debug >= 1) write(6,*) 'logical: ', lret
      endif
      end subroutine parse_get_logical
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to store characters corresponding to key.
!! \param[in] key  Key
!! \param[out] str  Characters to be stored
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_by_raw_key(key, str, lerr)
      implicit none
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: str
      logical, intent(out) :: lerr
      logical :: lerrsta, lerrend
      integer(8) :: ksta, kend
      integer(4) :: i
      lerr = .false.
      lasterrmsg = ''
      str = ''
      call hash_lookup(vsta, trim(key), ksta, lerrsta)
      call hash_lookup(vend, trim(key), kend, lerrend)
      if (lerrsta .or. lerrend) then
        write(lasterrmsg, '(a,a)') 'Cannot find key: ', trim(key)
        lerr = .true.
        return
      endif
      call hash_delete(vsta, trim(key), lerrsta)
      call hash_delete(vend, trim(key), lerrend)
      if (lerrsta .or. lerrend) then
        write(lasterrmsg, '(a,a)') 'Cannot delete key: ', trim(key)
        lerr = .true.
        return
      endif
      if (len(str) < kend - ksta + 1) then
        lasterrmsg = 'Too short string'
        lerr = .true.
        return
      endif
      do i = 1, kend-ksta+1
        str(i:i) = allstr(i+ksta-1)
      enddo
      end subroutine parse_get_by_raw_key
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to specify number of presence when the same
!!    keywords exist.
!! \param[in] key  Key to be searched.
!! \param[in] n  Variable of specified number of presence.
!<
      subroutine parse_specify_key(key, n)
      implicit none
      character(len=*), intent(in) :: key
      integer(4), intent(in) :: n
      increment = .false.
      call parse_complement_keyword(key, n)
      lastkey = trim(lastkey) // '/'
      end subroutine parse_specify_key
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to set keyword to be ignored.
!! \param[in] key  Key to be ignored.
!<
      subroutine parse_ignore_key(key)
      implicit none
      character(len=*), intent(in) :: key
      integer(8):: one=1
      logical:: lerr
      if (.not. use_ignored_key) then
        call hash_new(ignored_key)
        use_ignored_key = .true.
      endif
      call hash_insert(ignored_key, key, one, lerr)
      end subroutine parse_ignore_key
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to set characters corresponding to key.
!! \param[in] key  Key to be searched.
!! \param[out] str  Characters to be set.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_get_string(key, str, lerr)
      implicit none
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: str
      logical, intent(out) :: lerr
      call parse_complement_keyword(key)
      increment = .true.
      call parse_get_by_raw_key(trim(lastkey), str, lerr)
      if (debug >= 1) then
        if (lerr) then
          write(6,*) 'NG: ', trim(lastkey)
        else
          write(6,*) 'OK: ', trim(lastkey), ' ', trim(str)
        endif
      endif
      lastkey = trim(lastkey) // '/'
      end subroutine parse_get_string
! ---------------------------------------------------------------------
!>
!! \brief  (Subroutine to convert key given by user to internal form)
!! \param[in] origkey  Key given by user.
!! \param[in] n  Number of order of key.
!<
      subroutine parse_complement_keyword(origkey, n)
      implicit none
      character(len=*), intent(in) :: origkey
      integer(4), intent(in), optional :: n
      character(len=maxkey) :: key
      character(len=8) :: snum
      integer(4) :: nsl, nsta, nend, msl, msta, mend, mlen, minsl
      integer(4) :: nth, i, j
!     An example of 'origkey' is '/abc/cde/fghi'.
!     An example of 'key' is '/abc/cde/fghi/'.
!     An example of 'lastkey' is '/abc*2/cde*1/fghi*10/'.
!                                 |  +-mend           +-mlen
!                                 +--msta
!                                 msl : number of slashes
      lasterrmsg = ''
      key = trim(origkey) // '/'
      nsl = nslash(key)
      msl = nslash(lastkey)
      minsl = min(nsl, msl)
      nsta = 1
      msta = 1
      mlen = len_trim(lastkey)
      do i=2, minsl
        nend = nsta + scan(key(nsta+1:), '/') - 1
        mend = msta + scan(lastkey(msta+1:mlen), '*') - 1
        if (key(nsta+1:nend) /= lastkey(msta+1:mend)) then
          lastkey = lastkey(1:msta)
          msl = nslash(lastkey)
          minsl = min(nsl, msl)
          mlen = msta
          exit
        endif
        nsta = nend + 1
        msta = msta + scan(lastkey(msta+1:mlen), '/')
      enddo
      if (msta /= mlen) then
        lastkey = lastkey(1:msta)
        msl = nslash(lastkey)
        minsl = min(nsl, msl)
        mlen = msta
      endif
      if (i-1 == nsl) then
        msta = scan(lastkey, '*', .true.) + 1
        if (increment) then
          read(lastkey(msta:mlen-1), '(i8)') nth
          call string2int(lastkey(msta:mlen-1), nth)
          call int2string(nth + 1, snum)
        else
          snum = lastkey(msta:mlen-1)
        endif
        if (present(n)) then
          if (n > 0) then
            call int2string(n, snum)
          endif
        endif
        lastkey(msta:) = snum
      else
        snum = '*1'
        if (present(n)) then
          if (n > 0) then
            call int2string(n, snum)
            snum = '*' // snum
          endif
        endif
        do j=i, nsl-1
          nend = nsta + scan(key(nsta+1:), '/') - 1
          lastkey = lastkey(:mlen) // key(nsta+1:nend) // '*1/'
          mlen = mlen + nend - nsta + 3
          nsta = nend + 1
        enddo
        nend = nsta + scan(key(nsta+1:), '/') - 1
        lastkey = lastkey(:mlen) // key(nsta+1:nend) // snum
      endif
    
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (Internal subroutine to count slash.)
!! \param[in] string  Characters to be count.
!! \return  Number of slash.
!<
      integer(4) function nslash(string)
      implicit none
      character(len=*) :: string
      character:: c
      integer(4) :: i
      nslash = 0
      do i=1, len(string)
        c = key(i:i)
        if (c == '/') then
          nslash = nslash + 1
        else if (ichar(c) < 32) then
!         Only control character, not space(32).
          exit
        endif
      enddo
      end function nslash
      end subroutine parse_complement_keyword
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to finalize.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_close(lerr)
      use array
      implicit none
      logical, intent(out) :: lerr
      integer(4) :: iunused, i
      character(len=maxkey) :: key
      type(array_string) :: keys
      if (lallstr > 0) then
        deallocate(allstr)
        lallstr = -1
      endif
      call array_new(keys)
      call hash_keys(vsta, keys, lerr)
      iunused = array_number(keys)
      if (iunused > 0) then
        lerr = .true.
        lasterrmsg = 'There are unused items:'
        do i=1, iunused
          call array_pop(keys, key)
          lasterrmsg = trim(lasterrmsg) // ' ' // key
        enddo
      else
        lerr = .false.
      endif
      call array_destroy(keys)
      call hash_destroy(vsta)
      call hash_destroy(vend)
      call array_destroy(rawkeywords)
      if (use_ignored_key) then
        call hash_destroy(ignored_key)
        use_ignored_key=.false.
      endif
      end subroutine parse_close
! ---------------------------------------------------------------------
!>
!! \brief  Function to show detailed error message about parse.
!! \return  Characters of error message.
!<
      character(len=256) function parse_show_error_message()
      implicit none
      parse_show_error_message = lasterrmsg
      return
      end function parse_show_error_message
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to set debug mode.
!! \param[in] d  Logical value of debug mode.
!<
      subroutine parse_set_debug_logical(d)
      logical, intent(in) :: d
      if (d) then
        debug = 1
      else
        debug = 0
      endif
      end subroutine parse_set_debug_logical
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to set debug level.
!! \param[in] level  Integer value of debug level.
!<
      subroutine parse_set_debug_integer4(level)
      integer(4), intent(in) :: level
      debug = level
      end subroutine parse_set_debug_integer4
! =====================================================================

#ifdef GROEXT
! ---------------------------------------------------------------------
!>
!! \brief  Subroutine to read and parse input string.
!! \param[in] str Input string.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_open_string(str, lerr)
      use mpi_tool
      implicit none
      character(LEN=*), intent(in) :: str
      logical, intent(out) :: lerr
      integer(4) :: istatus
      character(len=128) :: s
      character(len=9):: form_type
      integer :: ionum_last
#ifdef MPI
      include 'mpif.h'
      integer(4) :: ierr, myrank
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
#endif  /* MPI */
      lastkey = '/'
      lallstr = -1
      lerr = .false.
      increment=.false.
      mode_ignored=.false.
      call array_new(rawkeywords)
      call array_push(rawkeywords, '')
#ifdef MPI
      if (myrank == 0) then
#endif  /* MPI */
        lasterrmsg = ''
        ionum_last = ionum
        call parse_open_string__load(str)
        call array_pop(rawkeywords, s)
        if (array_number(rawkeywords) > 0) then
          call array_join('/', rawkeywords, s)
          write(lasterrmsg, '(a)') 'Tag is not closed.  ' // s
          lerr = .true.
          return
        endif
#ifdef MPI
      endif
      call MPI_BCAST(lallstr,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(allstr, lallstr, MPI_CHARACTER, &
     &               0,MPI_COMM_WORLD,ierr)
#endif  /* MPI */
      call parse_analyze(lerr)
      contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (Internal subroutine) Count number of valid letters in the
!!    specified unit.
!<
      subroutine parse_open_string__load(str)
      use array
      implicit none
      character(LEN=*), intent(in) :: str
#ifdef __GFORTRAN__
      integer(4), parameter:: size_line=1024
#else /* __GFORTRAN__ */
      integer(4), parameter:: size_line=1024
#endif /* __GFORTRAN__ */
      integer(4) :: istatus, si, so, i, n_output, ic, nc
      integer(8):: lsize,lsizenew
      character, pointer :: allstrnew(:)
      character :: c
      character(len=128):: output
      character(len=size_line):: line
      logical :: ignore
      istatus = 0
      si = 0
      lallstr = 0
      ignore = .false.
      lsize = 10000
      allocate(allstr(lsize))
      do ic = 1, len(str)
          c = str(ic:ic)
          if (ichar(c) == 10) then
            ignore = .false.
          endif
          if (ignore) then
            cycle
          endif
          if (ichar(c) <= 32) then
            c = ' '
          endif
          if (c == '#') then
            ignore = .true.
            cycle
          else if (c == ';' .or. c == ',') then
            c = ' '
          endif
          call parse_open__state_value(c, si, so, output, n_output)
          if (lallstr + n_output > lsize) then
            lsizenew = lsize + lsize
            do while (lallstr + n_output > lsizenew)
              lsizenew = lsizenew + lsizenew
            enddo
            allocate(allstrnew(lsizenew))
            allstrnew(:lsize) = allstr(:lsize)
            deallocate(allstr)
            allstr => allstrnew
            lsize = lsizenew
          endif
          do i = 1, n_output
            allstr(lallstr+i) = output(i:i)
          enddo
          lallstr = lallstr + n_output
          si = so
      enddo
      if (debug >= 2) then
        call parse_write_allstr(allstr, lallstr)
      endif
      end subroutine parse_open_string__load
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (internal subroutine) Return an input state.
!! \param[in] c  Current letter to parse.
!! \param[in] si  Last state before parse.
!! \param[out] so  Result of state after parse.
!! \param[out] output  Decided string.
!! \param[out] n_output Length of decided string.
!! \details
!!   Table of state is following.
!!   \code
!!
!!  Consider 16 patterns in inputs.
!!          <a>  </a> val  k=v
!!   <a>     o     o    o    o
!!   </a>    o     o    o    o
!!   val     o     o    o    o
!!   k=v     o     o    o    o
!!
!!                     so
!!          0 10 11 12 13 20 21 22 23 99
!!      0  ss kv -- -- -- k- -- -- --  o    (Empty)
!!     10  -- ss kv K- -- kV -- -- --  o    a
!!     11  -- -V ss K- -- kV -- -- --  o    a (space)
!!     12  -- -- -- ss -v -- -- -- --  o    a=
!! si  13  -V -- -- -- ss kV -- -- --  o    a=A
!!     20  -- -- -- -- -- -- k- k- --  o    <
!!     21  -K -- -- -- -- -- ss -- --  o    <a
!!     22  -- -- -- -- -- -- -- -- k-  o    </
!!     23  -K -- -- -- -- -- -- -- ss  o    </a
!!
!!   k: Not decided keyword
!!   K: Decided keyword
!!   v: Not decided value
!!   V: Decided value
!!  ss: Recurrence
!!   \endcode
!<
      subroutine parse_open__state_value(c, si, so, output, n_output)
      use array
      implicit none
      character, intent(in)::c
      integer(4), intent(in):: si
      integer(4), intent(out):: so
      character(len=128), intent(out):: output
      integer(4), intent(out):: n_output
      character(len=72),save:: str=''
      character(len=128):: s
      integer(4),save:: length_str=0
      logical:: keep, lerr
!
      keep = .false.
      output = ''
      n_output = 0
      if (si == 0) then
        if (c == ' ') then
          so = 0
          return
        else if (normal_letter(ichar(c))) then
          so = 10
          keep = .true.
        else if (c == '<') then
          so = 20
        else
          so = 99
        endif
      else if (si == 10) then
        if (c == ' ') then
          so = 11
          keep = .true.
        else if (c == '=') then
          so = 12
!         Decide keyword
          call array_push(rawkeywords, str)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          if (debug >= 1) then
            write(6,*) trim(s), mode_ignored
          endif
          if (mode_ignored) then
            output = ''
            n_output = 0
            keep = .true.
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (c == '<') then
          so = 20
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (normal_letter(ichar(c))) then
          so = 10
          keep = .true.
        else
          so = 99
        endif
      else if (si == 11) then
        if (c == ' ') then
          so = 11
          return
        else if (c == '=') then
          so = 12
!         -- Decide keyword --
          call array_push(rawkeywords, str)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          if (debug >= 1) then
            write(6,*) trim(s), mode_ignored
          endif
          if (mode_ignored) then
            output = ''
            n_output = 0
            keep = .true.
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (c == '<') then
          so = 20
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
        else if (normal_letter(ichar(c))) then
          so = 10
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
          keep = .true.
        else
          so = 99
        endif
      else if (si == 12) then
        if (c == ' ') then
          so = 12
          return
        else if (normal_letter(ichar(c))) then
          so = 13
          keep = .true.
        else
          so = 99
        endif
      else if (si == 13) then
        if (c == ' ') then
          so = 0
!         -- Decide value --
          if (debug >= 1)  write(6,*) 'V: "', str(:length_str), '"'
          if (mode_ignored) then
            output = ''
            n_output = 0
          else
            output = ' ' // str
            n_output = 1 + length_str
          endif
          str = ''
          length_str = 0
          call array_pop(rawkeywords, s)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
        else if (normal_letter(ichar(c))) then
          so = 13
          keep = .true.
        else
          so = 99
        endif
      else if (si == 20) then
        if (c == '/') then
          so = 22
        else
          so = 21
          keep = .true.
        endif
      else if (si == 21) then
        if (c == '>') then
          so = 0
!         -- Decide keyword --
          call array_push(rawkeywords, str)
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          if (debug >= 1) then
            write(6,*) trim(s), mode_ignored
          endif
          output = str
          n_output = length_str
          str = ''
          length_str = 0
        else
          so = 21
          keep = .true.
        endif
      else if (si == 22) then
        if (c == '>') then
          so = 99
        else
          so = 23
          keep = .true.
        endif
      else if (si == 23) then
        if (c == '>') then
          so = 0
!         -- Decide keyword --
          call array_pop(rawkeywords, s)
          if (s(:length_str) .ne. str(:length_str)) then
            write(6,*) 'Warning: Keyword differ.'
            write(6,*) ' Start: "' // trim(s) // '".'
            write(6,*) ' End: "' // str(:length_str) // '".'
          endif
          call array_join('/', rawkeywords, s)
          if (use_ignored_key) then
            mode_ignored = hash_defined(ignored_key, s, lerr)
          endif
          output = str
          n_output = length_str
          str = ''
          length_str = 0
        else
          so = 23
          keep = .true.
        endif
      endif
      if (keep) then
        length_str = length_str + 1
        str(length_str:length_str) = c
      else
        n_output = n_output + 1
        output(n_output:n_output) = c
      endif
      if (debug >= 2) then
        write(6,*) c, si, so, '"', str(:length_str), '" "', &
     &             output(:n_output), '"'
      endif
      end subroutine parse_open__state_value
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>
!! \brief  (Intenal subroutine) Parse tags.
!! \param[out] lerr  True if error, otherwise False.
!<
      subroutine parse_analyze(lerr)
      use array
      use hash
      implicit none
      logical, intent(out) :: lerr
      character(len=128) :: s,keyword,merged,lastkeyword
      character(len=8) :: snum
      logical :: flag_equal=.false.
      type(array_string) :: indexedkeywords
      type(hash_i8) :: keyvarray
      integer(8) :: m, msta, mend, ksta, kend, is, ie, p, peq
      lerr = .true.
      call array_new(indexedkeywords)
      call array_push(indexedkeywords, '')
      call hash_new(keyvarray)
      call hash_new(vsta)
      call hash_new(vend)
      is = 1
      ie = lallstr
      do while (is > 0)
        is = parse_verify(allstr, is, ie, ' ')
        if (is == 0) then
!         case : End with space.
          exit
        endif
        if (allstr(is) == '<') then
          if (is+1 > ie) then
            lasterrmsg = 'End with <.'
            goto 900
          endif
          msta = is
          mend = parse_scan(allstr, is, ie, '>')
          if (mend == 0) then
            lasterrmsg = 'Parenthesis error.'
            goto 900
          endif
          if (allstr(msta+1) == '/') then
!           case : </keyword>
            ksta = msta + 2
            kend = mend - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            if (flag_equal) then
              call array_pop(indexedkeywords, s)
              call array_pop(rawkeywords, s)
              flag_equal = .false.
            endif
            call array_pop(indexedkeywords, s)
            call array_pop(rawkeywords, s)
            if (trim(s) /= trim(keyword)) then
              write(lasterrmsg, '(a,a,a,a)') &
     &            'Parenthesis keyword differ.', &
     &            trim(s), ' ', trim(keyword)
              goto 900
            endif
          else
!           case : <keyword>
            ksta = msta + 1
            kend = mend - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            if (flag_equal) then
              call array_pop(indexedkeywords, s)
              call array_pop(rawkeywords, s)
              flag_equal = .false.
            endif
            call array_join('/', indexedkeywords, merged)
            if (ichar(merged(2:2)) < 32) then
              merged = '/' // keyword
            else
              merged = trim(merged) // '/' // keyword
            endif
            call hash_lookup(keyvarray, merged, m, lerr)
            m = m + 1
            call hash_insert(keyvarray, merged, m, lerr)
            call int2string(m, snum)
            s = trim(keyword) // '*' // snum
            call array_push(indexedkeywords, s)
            call array_push(rawkeywords, keyword)
            lastkeyword = keyword
          endif
          is = mend + 1
        else
          p = parse_scan(allstr, is, ie, '<>= ')
          if (allstr(p) == ' ') then
!           Remove extra spaces before equal
            peq = parse_verify(allstr, p, ie, ' ')
            if (allstr(peq) == '=') then
              p = peq
            endif
          endif
          if (p == 0) then
!           case : End with value string.
            ksta = is
            kend = ie
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            is = 0
          else if (allstr(p) == '>') then
            lasterrmsg = 'Value contains >.'
            goto 900
          else if (allstr(p) == '<' .or. allstr(p) == ' ') then
!           case : value
            ksta = is
            kend = p - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            is = p
            call array_pop(indexedkeywords, s)
            call array_pop(rawkeywords, s)
            call array_join('/', indexedkeywords, merged)
            if (ichar(merged(2:2)) < 32) then
              merged = '/' // lastkeyword
            else
              merged = trim(merged) // '/' // lastkeyword
            endif
            call hash_lookup(keyvarray, merged, m, lerr)
            call hash_insert(keyvarray, merged, m+1, lerr)
            call int2string(m, snum)
            s = trim(lastkeyword) // '*' // snum
            call array_push(indexedkeywords, s)
            call array_push(rawkeywords, lastkeyword)
            call array_join('/', indexedkeywords, s)
            call hash_insert(vsta, s, ksta, lerr)
            call hash_insert(vend, s, kend, lerr)
          else if (allstr(p) == '=') then
!           case : keyword=
            ksta = is
            kend = p - 1
            call chars2string(ksta, kend, keyword, lerr)
            if (lerr) goto 900
            is = p + 1
            if (flag_equal) then
              call array_pop(indexedkeywords, s)
              call array_pop(rawkeywords, s)
              flag_equal = .false.
            endif
            call array_join('/', indexedkeywords, merged)
            if (ichar(merged(2:2)) < 32) then
              merged = '/' // keyword
            else
              merged = trim(merged) // '/' // keyword
            endif
            call hash_lookup(keyvarray, merged, m, lerr)
            m = m + 1
            call hash_insert(keyvarray, merged, m, lerr)
            call int2string(m, snum)
            s = trim(keyword) // '*' // snum
            call array_push(indexedkeywords, s)
            call array_push(rawkeywords, keyword)
            lastkeyword = keyword
            flag_equal = .true.
          endif
        endif
      enddo
      call array_pop(indexedkeywords, s)
      if (array_number(indexedkeywords) /= 0) then
        lasterrmsg = 'ERROR: Unexpected EOF.'
        goto 900
      endif
      lerr = .false.
  900 continue
      call array_destroy(indexedkeywords)
      call array_destroy(rawkeywords)
      call hash_destroy(keyvarray)
      return
      end subroutine parse_analyze
      end subroutine parse_open_string

#endif /* GROEXT */

      end module parse
