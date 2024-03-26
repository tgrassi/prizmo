module prizmo_utils
  use prizmo_commons
contains

  ! ************************
  function get_electrons(x) result(ne)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::ne

    !! PREPROCESS_ELECTRONS
    !! PREPROCESS_END

    if(ne < 0d0) then
      print *, "ERROR: negative electrons!"
      stop
    end if

    ne = max(ne, 0d0)

  end function get_electrons

  ! ************************
  function get_rho(x) result(rho)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::rho

    !! PREPROCESS_GET_RHO
    !! PREPROCESS_END

  end function get_rho

  ! ************************
  function get_Hnuclei(x) result(nH)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nH

    !! PREPROCESS_HNUCLEI
    !! PREPROCESS_END

    if(nH < 0d0) then
      print *, "ERROR: negative H nuclei!"
      print *, nH
      !stop
    end if

    nH = max(nH, 0d0)

  end function get_Hnuclei

  ! ************************
  function get_Cnuclei(x) result(nC)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nC

    !! PREPROCESS_CNUCLEI
    !! PREPROCESS_END

    if(nC < 0d0) then
      print *, "ERROR: negative C nuclei!"
      print *, nC
      !stop
    end if

    nC = max(nC, 0d0)

  end function get_Cnuclei

  ! ************************
  function get_Onuclei(x) result(nO)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nO

    !! PREPROCESS_ONUCLEI
    !! PREPROCESS_END

    if(nO < 0d0) then
      print *, "ERROR: negative O nuclei!"
      print *, nO
      !stop
    end if

    nO = max(nO, 0d0)

  end function get_Onuclei

  ! ************************
  function get_HEnuclei(x) result(nHE)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nHE

    !! PREPROCESS_HENUCLEI
    !! PREPROCESS_END

    if(nHE < 0d0) then
      print *, "ERROR: negative HE nuclei!"
      print *, nHE
      !stop
    end if

    nHE = max(nHE, 0d0)

  end function get_HEnuclei

  ! **************************
  subroutine ranker(array, nbest, array_string)
    implicit none
    real*8,intent(in):: array(:)
    integer,intent(in)::nbest
    character(len=*),intent(in):: array_string(:)
    integer::idxs(size(array)), i
    real*8::amax

    idxs = argsort_r(array)

    amax = array(idxs(1))

    print *, "**********************"
    do i=1,min(nbest, size(array))
      print '(2I5,2E17.8e3,x,a30)', i, idxs(i), array(idxs(i)), array(idxs(i)) / (amax + 1d-40), trim(array_string(idxs(i)))
    end do

  end subroutine ranker

  ! **********************
  function argsort_r(array) result(idxs)
    implicit none
    real*8,intent(in)::array(:)
    integer::idxs(size(array)), na

    na = size(array)

    idxs = argsort(array)
    idxs = idxs(na:1:-1)

  end function

  ! **********************
  function argsort(array) result(idxs)
    implicit none
    real*8,intent(in)::array(:)
    real*8::a(size(array)), rtmp
    integer::idxs(size(array)), i, na, itmp
    logical::swap

    na = size(array)
    do i=1,na
      idxs(i) = i
    end do

    a = array

    do
      swap = .false.
      do i=2,na
        if(a(i-1) > a(i)) then
          rtmp = a(i)
          a(i) = a(i-1)
          a(i-1) = rtmp
          itmp = idxs(i)
          idxs(i) = idxs(i-1)
          idxs(i-1) = itmp
          swap = .true.
        end if
      end do
      if(.not.swap) exit
    end do

  end function argsort

end module prizmo_utils
