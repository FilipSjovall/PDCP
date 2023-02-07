!!!
module mesh
! Program that reads from ASCII-mesh file
!!!
!
!
! The file hej.msh is in ASCII - 1 format
!
!
!
!
contains 


subroutine read_mesh(fname, coord, enod)
    !character(len=100) :: fname ="~/Documents/Kurser/PDC Summer School/Project/hej.msh"
    character(len=*) :: fname
    integer :: n, iostat, nstart, nend, elstart, elend, idum, idum2, ierr
    integer, dimension(1,4) :: dummy_int
    real(kind=8) :: x
    character(len=100) :: dummy_char
    real(kind=8), allocatable, intent(inout) :: coord(:,:)
    integer, allocatable, intent(inout) :: enod(:,:)

    print * , '------------------ READING MESH ------------------------'
    n = 0
    open(45, file=fname, status="old",action="read")
    !
    !
    ! Loop to find the number of lines in the file
    do
        read(45 ,* ,IOSTAT=iostat) dummy_char
        n = n+1
        if (dummy_char(1:6) =="$NOD" ) then
            nstart = n
        elseif ( dummy_char(1:9)=="$ENDNOD" ) then
            nend   = n
        elseif ( dummy_char(1:9)=="$ELM" ) then
            elstart = n
        elseif ( dummy_char(1:13)=="$ENDELM" ) then
            elend   = n
        endif
        if(iostat.lt.0) then
            exit
        endif
    end do
    !
    print *, "Number of nodes ", nend-nstart-2

    print *, "Number of elements ", elend-elstart-3

    1 rewind(45)


    if(allocated(coord)) deallocate(coord)
    if(allocated(enod)) deallocate(enod)

    allocate(coord(nend-nstart-2,2),stat=ierr)

    allocate(enod(elend-elstart-3,7),stat=ierr)
    ! Loop to find the number of lines in the file
    do n=1,elend
        ! Sista villkoret borde vara första if-statementet
        if (( n>nstart+1).and.(n.le.nend-1)) then
            read(45 ,* ,IOSTAT=iostat), idum2, coord(n-nstart-1,1), coord(n-nstart-1,2), idum
            !print *, n-nstart-1, coord(n-nstart-1,:)
        elseif ((n > elstart+1).and.(n<elend-1)) then
            read(45 ,* ,IOSTAT=iostat), enod(n-elstart-1,1), dummy_int, enod(n-elstart-1,2:7)
            !print *, enod(n-elstart-1,:), "\\"
        else
            read(45 ,* ,IOSTAT=iostat) dummy_char
            !print *, dummy_char
        endif
    end do
end subroutine read_mesh

end module mesh