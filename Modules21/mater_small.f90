module mater_small

! last modified
! M. Ristinmaa 2011-04-28
!  - module changed name
! M. Ristinmaa 2011-05-02
!  - introduced 4x4 matrix plane strain hooke
!  - mises implemented
! ------------------------------------------------------------------------------
!

implicit none

interface mises
   module procedure mises1
   module procedure mises2
end interface


!-------------------------------------------------------------------------------

contains

subroutine mises1(svm,ss)
  implicit none
  double precision                :: svm, ss(:)
  double precision                :: s(3,3)

  if (size(ss).eq.3) then
    s(1,:)=(/ss(1),ss(3),0d0/)
    s(2,:)=(/ss(3),ss(2),0d0/)
    s(3,:)=(/0d0,0d0,0d0/)
  elseif (size(ss).eq.4) then
    s(1,:)=(/ss(1),ss(4),0d0/)
    s(2,:)=(/ss(4),ss(2),0d0/)
    s(3,:)=(/0d0,0d0,ss(3)/)
  else
   stop 'mises not implemented'
  endif

  svm=dsqrt(s(1,1)**2D0+s(2,2)**2D0+s(3,3)**2D0-s(1,1)*s(2,2)-s(1,1)*s(3,3)-s(2,2)*s(3,3)&
                         +3D0*s(1,2)**2D0+3D0*s(1,3)**2D0+3D0*s(2,3)**2D0)

  return
end subroutine mises1


subroutine mises2(svm,ss)
  implicit none
  double precision                :: svm(:), ss(:,:)
  double precision                :: s(3,3)
  integer                         :: ii

  do ii=1,size(ss,2)
    if (size(ss,1).eq.3) then
      s(1,:)=(/ss(ii,1),ss(ii,3),0d0/)
      s(2,:)=(/ss(ii,3),ss(ii,2),0d0/)
      s(3,:)=(/0d0,0d0,0d0/)
    else
      stop 'mises not implemented'
    endif

    svm(ii)=dsqrt(s(1,1)**2D0+s(2,2)**2D0+s(3,3)**2D0-s(1,1)*s(2,2)-s(1,1)*s(3,3)-s(2,2)*s(3,3)&
                         +3D0*s(1,2)**2D0+3D0*s(1,3)**2D0+3D0*s(2,3)**2D0)
  end do

  return
end subroutine mises2


subroutine hooke(D,ptype,E,v)
  implicit none
  double precision        :: D(:,:), E, v
  integer                 :: ptype

  if (ptype.eq.1) then
        D=E/(1d0-v**2d0)*reshape((/1d0,   v, 0d0,&
                                   v  , 1d0, 0d0,&
                                   0d0, 0d0, (1d0-v)*0.5d0/),&
                                  [3,3], order=[2,1])
  elseif (ptype.eq.2) then
        D=E/(1d0+v)/(1d0-2d0*v)&
                       *reshape((/1d0-v,     v,   0d0,&
                                  v    , 1d0-v,   0d0,&
                                  0d0  ,   0d0,  (1d0-2d0*v)*0.5d0/), &
                                  [3,3], order=[2,1])
  else
    write(*,*) 'Error ! Check first argument, ptype=1 or 2 allowed'
    stop
  end if

return
end subroutine hooke

	
	
end module mater_small
