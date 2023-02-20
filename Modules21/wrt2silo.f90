module wrt2silo

! =======================================================================================================
!
!   Module wrt2silo contains subroutines for writing Fortran data to Silo files.
!   The routines require the Silo libraries to be installed. These can be obtained from:
!   https://wci.llnl.gov/codes/visit/3rd_party
!
!   Included routines:
!
!      openSilo                   Open a Silo file and assign a pointer and a filename
!      closeSilo                  Close a Silo file
!      mkQmeshCollinear           Creates a quad mesh collinear with the coordinate axes (2D or 3D)
!      mkQmesh                    Creates a quad mesh collinear or non-collinear with the coordinate axes
!                                 (2D or 3D). Similar as 'mkQmeshCollinear', but based on node coordinate
!                                 input.
!      plotSolidMesh              Creates an unstructured mesh (2D or 3D)
!      plotVar                    Write nodal or zonal values to an unstructured mesh (scalar integer
!                                 or double precision values are supported)
!      mkQmeshNonCollinear        Creates a quad mesh non-collinear with the coordinate axes (2D or 3D)
!      setQmeshZoneVar            Assigns a zone variable to a quad mesh (2D or 3D, integer or double)
!      plotGraph                  Creates a xy-curve (integer or double)
!      mk3dQuadMultiMeshCollinear Creates a 3D quad multimesh, collinear with the coordinate axes.
!
!      plotPointMesh              Creates a 2D or 3D point mesh
!
!   Last rev.:
!      H. Hallberg, 2011-05-11
!
! =======================================================================================================

use silowrap

implicit none

!include "silo.inc"

! Make module methods inaccessible to outside calls
! =================================================
private :: setQmeshZoneDVar3d, setQmeshZoneIVar3d
private :: setQmeshZoneDVar2d, setQmeshZoneIVar2d
private :: mkQmesh2dCollinear, mkQmesh3dCollinear
private :: mkQmesh2dNonCollinear, mkQmesh3dNonCollinear
private :: mkQmesh2dNonCollinearStrain, mkQmesh3dNonCollinearStrain
private :: setDCurve1, setICurve1, setDCurve2, setICurve2
private :: mkQmesh3dShiftedCollinear
private :: setUmeshDVar, setUmeshIVar, setUmeshDVar2, setUmeshIVar2
private :: createSiloDb1, createSiloDb2
private :: closeSiloDb1, closeSiloDb2

! Declare accessible module interfaces
! ====================================
interface setQmeshZoneVar
   module procedure setQmeshZoneDVar3d
   module procedure setQmeshZoneIVar3d
   module procedure setQmeshZoneDVar2d
   module procedure setQmeshZoneIVar2d
end interface

interface mkQmeshCollinear
   module procedure mkQmesh3dCollinear
   module procedure mkQmesh2dCollinear
   module procedure mkQmesh3dShiftedCollinear
   module procedure mkQmesh3dShiftedCollinearStrain
end interface

interface plotVar
   module procedure setUmeshDVar
   module procedure setUmeshIVar
   module procedure setUmeshDVar2
   module procedure setUmeshIVar2
end interface

interface mkQmeshNonCollinear
   module procedure mkQmesh2dNonCollinear
   module procedure mkQmesh2dNonCollinearStrain
   module procedure mkQmesh3dNonCollinear
   module procedure mkQmesh3dNonCollinearStrain
end interface

interface plotGraph
   module procedure setDCurve1
   module procedure setDCurve2
   module procedure setICurve1
   module procedure setICurve2
end interface

interface openSilo
   module procedure createSiloDb1
   module procedure createSiloDb2
end interface

interface closeSilo
   module procedure closeSiloDb1
   module procedure closeSiloDb2
end interface

contains

! =======================================================================================================

subroutine createSiloDb1(dbfile,siloname)
   ! ===========================================================
   ! This routine opens Silo file.
   !
   !    Input:
   !          siloname          Name of the new Silo file
   !
   !    Output:
   !          dbfile            Pointer to the open Silo file
   !          
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================
   
   implicit none

   character(*), intent(in) :: siloname
   integer, intent(out)     :: dbfile
   integer                  :: ierr

   ! Create Silo file
   ! ================
   ierr = dbcreate(trim(siloname), len_trim(siloname), db_clobber, db_local, db_f77null, 0, db_pdb, dbfile)
   if (dbfile.eq.-1) then
      write(*,*) 'Could not create Silo file ', trim(siloname),'. Program execution stopped!'
      stop
   end if

   return
   
end subroutine createSiloDb1

! =======================================================================================================

subroutine createSiloDb2(dbfile,outfileid,siloname)
   ! ===========================================================
   ! This routine opens Silo file.
   !
   !    Input:
   !          outfileid         File-ID for program output
   !          siloname          Name of the new Silo file
   !
   !    Output:
   !          dbfile            Pointer to the open Silo file
   !          
   ! Last rev.:
   !    2011-04-07: H. Hallberg
   ! ===========================================================
   
   implicit none

   integer, intent(in)      :: outfileid
   character(*), intent(in) :: siloname
   integer, intent(out)     :: dbfile
   integer                  :: ierr

   ! Create Silo file
   ! ================
   ierr = dbcreate(trim(siloname), len_trim(siloname), db_clobber, db_local, db_f77null, 0, db_pdb, dbfile)
   if (dbfile.eq.-1) then
      write(outfileid,*) 'Could not create Silo file ', trim(siloname),'. Program execution stopped!'
      stop
   end if

   return
   
end subroutine createSiloDb2

! =======================================================================================================

subroutine closeSiloDb1(dbfile)
   ! ===========================================================
   ! This routine closes an open Silo file.
   !
   !    Input:
   !          dbfile            Pointer to the open Silo file
   !
   !    Output:
   !          
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================
   
   implicit none

   integer, intent(in) :: dbfile
   integer             :: ierr

   ierr = dbclose(dbfile)
   if (ierr.eq.-1) then
      write(*,*) 'Could not close Silo file. Program execution stopped!'
      stop
   end if

   return
   
end subroutine closeSiloDb1

! =======================================================================================================

subroutine closeSiloDb2(dbfile,outfileid)
   ! ===========================================================
   ! This routine closes an open Silo file.
   !
   !    Input:
   !          dbfile            Pointer to the open Silo file
   !          outfileid         File-ID for program output 
   !
   !    Output:
   !          
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================
   
   implicit none

   integer, intent(in) :: dbfile, outfileid
   integer             :: ierr

   ierr = dbclose(dbfile)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not close Silo file. Program execution stopped!'
      stop
   end if

   return
   
end subroutine closeSiloDb2

! =======================================================================================================

subroutine mkQmesh3dCollinear(dbfile,meshname,outfileid,dim1,dim2,dim3,ndiv1,ndiv2,ndiv3)
   ! ===========================================================
   ! This routine creates a 3D quad mesh collinear with the
   ! coordinate axes.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          dim1,dim2,dim3    3D dimensions of the domain
   !          ndiv1,ndiv2,ndiv3 Number of element/cell
   !                            divisions along the three
   !                            spatial dimensions
   !
   !    Output:
   !          
   !          
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3
   double precision, intent(in) :: dim1, dim2, dim3
   character(*), intent(in)     :: meshname
   integer                      :: ierr, err, ndims, dims(3), i
   double precision             :: xc(ndiv1+1), yc(ndiv2+1), zc(ndiv3+1), dx, dy, dz

   ! Write quad mesh data
   ! ====================
   ndims = 3
   dims  = (/ ndiv1+1, ndiv2+1, ndiv3+1 /)
   xc = 0D0
   yc = 0D0
   zc = 0D0
   dx = dim1/ndiv1
   dy = dim2/ndiv2
   dz = dim3/ndiv3
   do i=1,ndiv1
      xc(i+1) = i*dx
   end do
   do i=1,ndiv2
      yc(i+1) = i*dy
   end do
   do i=1,ndiv3
      zc(i+1) = i*dz
   end do
   ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,zc,dims,ndims,db_double,db_collinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write collinear 3D quad mesh ', trim(meshname), ' to Silo file. Program execution stopped!'
      stop
   end if
   
   return

end subroutine mkQmesh3dCollinear

! =======================================================================================================

subroutine mkQmesh2dCollinear(dbfile,meshname,outfileid,dim1,dim2,ndiv1,ndiv2)
   ! ===========================================================
   ! This routine creates a 2D quad mesh collinear with the
   ! coordinate axes.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          dim1,dim2         2D dimensions of the domain
   !          ndiv1,ndiv2       Number of element/cell
   !                            divisions along the two
   !                            spatial dimensions
   !
   !    Output:
   !          
   !          
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2
   double precision, intent(in) :: dim1, dim2
   character(*), intent(in)     :: meshname
   integer                      :: ierr, err, ndims, dims(2), i
   double precision             :: xc(ndiv1+1), yc(ndiv2+1), dx, dy
   
   ! Write quad mesh data
   ! ====================
   ndims = 2
   dims  = (/ ndiv1+1, ndiv2+1 /)
   xc = 0D0
   yc = 0D0
   dx = dim1/ndiv1
   dy = dim2/ndiv2
   do i=1,ndiv1
      xc(i+1) = i*dx
   end do
   do i=1,ndiv2
      yc(i+1) = i*dy
   end do
   ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,db_f77null, &                
                  dims,ndims,db_double,db_collinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write collinear 2D quad mesh ', trim(meshname), ' to Silo file. Program execution stopped!'
      stop
   end if

   return

end subroutine mkQmesh2dCollinear

! =======================================================================================================

subroutine setQmeshZoneDVar3d(dbfile,meshname,outfileid,ndiv1,ndiv2,ndiv3,dvar,nvar,varname,varunit,dtime,istep)
   ! ===========================================================
   ! This routine writes a zone-centered variable of Fortran
   ! type "double" to an existing Silo file with a defined quad
   ! mesh.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2,ndiv3 Number of cell divisions
   !                            along the three spatial
   !                            dimensions
   !          dvar              Variable array (type double)
   !                            to be written to the Silo file
   !          nvar              Number of variable values
   !          varname           Character string containing the
   !                            name of the variable
   !          varunit           Character string containing the
   !                            unit of the variable
   !          dtime             Current total solution time
   !          istep             Current solution step
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3, istep, nvar
   double precision, intent(in) :: dtime, dvar(nvar)
   character(*), intent(in)     :: varname, varunit, meshname
   integer                      :: dims(3), ndims, err, ierr, i, j, k, index, optlistid
   real                         :: zonal(ndiv1,ndiv2,ndiv3)
   
   ndims = 3
   dims  = (/ ndiv1, ndiv2, ndiv3 /)

   index = 1
   do k=1,ndiv3
      do j=1,ndiv2
         do i=1,ndiv1
            zonal(i,j,k) = dvar(index)
            index = index + 1
         end do
      end do
   end do
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(3,optlistid)
   err = dbaddiopt(optlistid, dbopt_cycle, istep)
   err = dbadddopt(optlistid, dbopt_dtime, dtime)
   err = dbaddcopt(optlistid, dbopt_units, trim(varunit), len_trim(varunit))
   
   err = dbputqv1(dbfile,varname,len_trim(varname),meshname, &
                  len_trim(meshname),zonal,dims,ndims,db_f77null,0,db_float,db_zonecent,optlistid,ierr)
   if (err.eq.-1) then
      write(outfileid,*) 'Could not write zone variable ',trim(varname), ' of type double to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
           
   return

end subroutine setQmeshZoneDVar3d

! =======================================================================================================

subroutine setQmeshZoneIVar3d(dbfile,meshname,outfileid,ndiv1,ndiv2,ndiv3,ivar,nvar,varname,varunit,dtime,istep)
   ! ===========================================================
   ! This routine writes a zone-centered variable of Fortran
   ! type "integer" to an existing Silo file with a defined quad
   ! mesh.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2,ndiv3 Number of cell divisions
   !                            along the three spatial
   !                            dimensions
   !          ivar              Variable array (type integer)
   !                            to be written to the Silo file
   !          nvar              Number of variable values
   !          varname           Character string containing the
   !                            name of the variable
   !          varunit           Character string containing the
   !                            unit of the variable
   !          dtime             Current total solution time
   !          istep             Current solution step
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2010-12-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3, istep, nvar
   integer, intent(in)          :: ivar(nvar)
   double precision, intent(in) :: dtime
   character(*), intent(in)     :: varname, varunit, meshname
   integer                      :: dims(3), ndims, err, ierr, i, j, k, index, optlistid
   integer                      :: zonal(ndiv1,ndiv2,ndiv3)

   ndims = 3
   dims  = (/ ndiv1, ndiv2, ndiv3 /)
   
   index = 1
   do k=1,ndiv3
      do j=1,ndiv2
         do i=1,ndiv1
            zonal(i,j,k) = ivar(index)
            index = index + 1
         end do
      end do
   end do
  
   ! Create option list
   ! ==================
   err = dbmkoptlist(3,optlistid)
   err = dbaddiopt(optlistid, dbopt_cycle, istep)
   err = dbadddopt(optlistid, dbopt_dtime, dtime)
   err = dbaddcopt(optlistid, dbopt_units, trim(varunit), len_trim(varunit))
  
   err = dbputqv1(dbfile,varname,len_trim(varname),meshname, &
         len_trim(meshname),zonal,dims,ndims,db_f77null,0,db_int,db_zonecent,optlistid,ierr)
   if (err.eq.-1) then
      write(outfileid,*) 'Could not write zone variable ', trim(varname), &
                         ' of type integer to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
   
   return

end subroutine setQmeshZoneIVar3d

! =======================================================================================================

subroutine setQmeshZoneDVar2d(dbfile,meshname,outfileid,ndiv1,ndiv2,dvar,nvar,varname,varunit,dtime,istep)
   ! ===========================================================
   ! This routine writes a zone-centered variable of Fortran
   ! type "double" to an existing Silo file with a defined quad
   ! mesh.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2       Number of cell divisions
   !                            along the two spatial
   !                            dimensions
   !          dvar              Variable array (type double)
   !                            to be written to the Silo file
   !          nvar              Number of variable values
   !          varname           Character string containing the
   !                            name of the variable
   !          varunit           Character string containing the
   !                            unit of the variable
   !          dtime             Current total solution time
   !          istep             Current solution step
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, istep, nvar
   double precision, intent(in) :: dtime, dvar(nvar)
   character(*), intent(in)     :: varname, varunit, meshname
   integer                      :: dims(2), ndims, err, ierr, i, j, index, optlistid
   real                         :: zonal(ndiv1,ndiv2)
   
   ndims = 2
   dims  = (/ ndiv1, ndiv2 /)
   
   index = 1
   do j=1,ndiv2
      do i=1,ndiv1
         zonal(i,j) = dvar(index)
         index = index + 1
      end do
   end do
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(3,optlistid)
   err = dbaddiopt(optlistid, dbopt_cycle, istep)
   err = dbadddopt(optlistid, dbopt_dtime, dtime)
   err = dbaddcopt(optlistid, dbopt_units, trim(varunit), len_trim(varunit))
   
   err = dbputqv1(dbfile,varname,len_trim(varname),meshname, &
                  len_trim(meshname),zonal,dims,ndims,db_f77null,0,db_float,db_zonecent,optlistid,ierr)
   if (err.eq.-1) then
      write(outfileid,*) 'Could not write zone variable ', trim(varname), ' of type double to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
   
   return

end subroutine setQmeshZoneDVar2d

! =======================================================================================================

subroutine setQmeshZoneIVar2d(dbfile,meshname,outfileid,ndiv1,ndiv2,ivar,nvar,varname,varunit,dtime,istep)
   ! ===========================================================
   ! This routine writes a zone-centered variable of Fortran
   ! type "double" to an existing Silo file with a defined quad
   ! mesh.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2       Number of cell divisions
   !                            along the two spatial
   !                            dimensions
   !          ivar              Variable array (type integer)
   !                            to be written to the Silo file
   !          nvar              Number of variable values
   !          varname           Character string containing the
   !                            name of the variable
   !          varunit           Character string containing the
   !                            unit of the variable
   !          dtime             Current total solution time
   !          istep             Current solution step
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, istep, nvar
   integer, intent(in)          :: ivar(nvar)
   double precision, intent(in) :: dtime
   character(*), intent(in)     :: varname, varunit, meshname
   integer                      :: dims(2), ndims, err, ierr, i, j, index, optlistid
   integer                      :: zonal(ndiv1,ndiv2)
   
   ndims = 2
   dims  = (/ ndiv1, ndiv2 /)
   
   index = 1
   do j=1,ndiv2
      do i=1,ndiv1
         zonal(i,j) = ivar(index)
         index = index + 1
      end do
   end do
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(3,optlistid)
   err = dbaddiopt(optlistid, dbopt_cycle, istep)
   err = dbadddopt(optlistid, dbopt_dtime, dtime)
   err = dbaddcopt(optlistid, dbopt_units, trim(varunit), len_trim(varunit))
   
   err = dbputqv1(dbfile,varname,len_trim(varname),meshname, &
                  len_trim(meshname),zonal,dims,ndims,db_f77null,0,db_int,db_zonecent,optlistid,ierr)
   if (err.eq.-1) then
      write(outfileid,*) 'Could not write zone variable ', trim(varname), &
                         ' of type integer to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
           
   return

end subroutine setQmeshZoneIVar2d

! =======================================================================================================

subroutine setICurve1(dbfile,xvals,yvals,xlabel,ylabel,curvename)
   ! ===========================================================
   ! This routine writes integer x-y data to an existing Silo
   ! file.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          xvals             x-values of the data points
   !          yvals             y-values of the data points
   !          xlabel            x-axis label (character string)
   !          ylabel            y-axis label (character string)
   !          curvename         Curve title (character string)
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)      :: dbfile
   integer, intent(in)      :: xvals(:), yvals(:)
   character(*), intent(in) :: xlabel, ylabel, curvename
   integer                  :: optlistid, err, ierr, npoints
   
   npoints = size(xvals,1) ! Number of data points
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(2,optlistid) 
   err = dbaddcopt(optlistid, dbopt_xlabel, xlabel, len_trim(xlabel))
   err = dbaddcopt(optlistid, dbopt_ylabel, ylabel, len_trim(ylabel))
     
   err = dbputcurve(dbfile, trim(curvename), len_trim(curvename), xvals, yvals, db_int, npoints, optlistid, ierr)

   if (err.eq.-1) then
      write(*,*) 'Could not write x-y curve ', trim(curvename), &
                 ' of type integer to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
           
   return

end subroutine setICurve1

! =======================================================================================================

subroutine setICurve2(dbfile,outfileid,xvals,yvals,xlabel,ylabel,curvename)
   ! ===========================================================
   ! This routine writes integer x-y data to an existing Silo
   ! file.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          outfileid         File-ID for program output 
   !          xvals             x-values of the data points
   !          yvals             y-values of the data points
   !          xlabel            x-axis label (character string)
   !          ylabel            y-axis label (character string)
   !          curvename         Curve title (character string)
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)      :: dbfile, outfileid
   integer, intent(in)      :: xvals(:), yvals(:)
   character(*), intent(in) :: xlabel, ylabel, curvename
   integer                  :: optlistid, err, ierr, npoints
   
   npoints = size(xvals,1) ! Number of data points
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(2,optlistid) 
   err = dbaddcopt(optlistid, dbopt_xlabel, xlabel, len_trim(xlabel))
   err = dbaddcopt(optlistid, dbopt_ylabel, ylabel, len_trim(ylabel))
     
   err = dbputcurve(dbfile, trim(curvename), len_trim(curvename), xvals, yvals, db_int, npoints, optlistid, ierr)

   if (err.eq.-1) then
      write(outfileid,*) 'Could not write x-y curve ', trim(curvename), &
                         ' of type integer to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
           
   return

end subroutine setICurve2

! =======================================================================================================

subroutine setDCurve1(dbfile,xvals,yvals,xlabel,ylabel,curvename)
   ! ===========================================================
   ! This routine writes float x-y data to an existing Silo
   ! file.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          xvals             x-values of the data points
   !          yvals             y-values of the data points
   !          xlabel            x-axis label (character string)
   !          ylabel            y-axis label (character string)
   !          curvename         Curve title (character string)
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile
   double precision, intent(in) :: xvals(:), yvals(:)
   character(*), intent(in)     :: xlabel, ylabel, curvename
   integer                      :: optlistid, err, ierr, npoints
   
   npoints = size(xvals,1) ! Number of data points
  
   ! Create option list
   ! ==================
   err = dbmkoptlist(2,optlistid) 
   err = dbaddcopt(optlistid, dbopt_xlabel, xlabel, len_trim(xlabel))
   err = dbaddcopt(optlistid, dbopt_ylabel, ylabel, len_trim(ylabel))

   err = dbputcurve(dbfile, trim(curvename), len_trim(curvename), xvals, yvals, db_double, npoints, optlistid, ierr)

   if (err.eq.-1) then
      write(*,*) 'Could not write x-y curve ', trim(curvename), ' of type integer to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
           
   return

end subroutine setDCurve1

! =======================================================================================================

subroutine setDCurve2(dbfile,outfileid,xvals,yvals,xlabel,ylabel,curvename)
   ! ===========================================================
   ! This routine writes float x-y data to an existing Silo
   ! file.
   !
   !    Input:
   !          dbfile            Silo file identifier (integer)
   !          outfileid         File-ID for program output 
   !          xvals             x-values of the data points
   !          yvals             y-values of the data points
   !          xlabel            x-axis label (character string)
   !          ylabel            y-axis label (character string)
   !          curvename         Curve title (character string)
   !
   !    Output:
   !                            Data is written to the
   !                            existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid
   double precision, intent(in) :: xvals(:), yvals(:)
   character(*), intent(in)     :: xlabel, ylabel, curvename
   integer                      :: optlistid, err, ierr, npoints
   
   npoints = size(xvals,1) ! Number of data points
  
   ! Create option list
   ! ==================
   err = dbmkoptlist(2,optlistid) 
   err = dbaddcopt(optlistid, dbopt_xlabel, xlabel, len_trim(xlabel))
   err = dbaddcopt(optlistid, dbopt_ylabel, ylabel, len_trim(ylabel))

   err = dbputcurve(dbfile, trim(curvename), len_trim(curvename), xvals, yvals, db_double, npoints, optlistid, ierr)

   if (err.eq.-1) then
      write(outfileid,*) 'Could not write x-y curve ', trim(curvename), ' of type integer to Silo file. Program execution stopped!'
      stop
   end if

   ! Free option list
   ! ================
   err = dbfreeoptlist(optlistid)
           
   return

end subroutine setDCurve2

! =======================================================================================================

subroutine mkQmesh2dNonCollinearStrain(dbfile,meshname,outfileid,ndiv1,ndiv2,xcoord,ycoord,epsx,epsy)
   ! ===========================================================
   ! This routine creates a 2D quad mesh in an existing Silo
   ! file. The mesh is not collinear with the coordinate axes.
   !
   !    Input:
   !          dbfile            ID of the created Silo file
   !                            (integer value)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2       Number of cell divisions
   !                            along the two spatial
   !                            dimensions
   !          xcoord            x-coordinates of the nodes
   !          ycoord            y-coordinates of the nodes
   !          epsx              Strain component in the x-dir.
   !          epsy              Strain component in the y-dir.
   !
   !    Output:
   !          
   !          
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2
   !double precision, intent(in) :: xcoord((ndiv1+1)*(ndiv2+1)), ycoord((ndiv1+1)*(ndiv2+1))
   double precision, intent(in) :: xcoord(:), ycoord(:), epsx, epsy
   character(*), intent(in)     :: meshname

   !double precision             :: xc(ndiv1+1,ndiv2+1), yc(ndiv1+1,ndiv2+1)
   real                         :: xc(ndiv1+1,ndiv2+1), yc(ndiv1+1,ndiv2+1)
   integer                      :: ierr, err, ndims, dims(2), dim1, dim2, i
 
   ! Write quad mesh data
   ! ====================
   ndims = 2
   dims  = (/ ndiv1+1, ndiv2+1 /)
   i     = 1
   do dim2=1,dims(2)
      do dim1=1,dims(1)
         xc(dim1,dim2) = xcoord(i)*epsx
         yc(dim1,dim2) = ycoord(i)*epsy
         i = i + 1
      end do
   end do
   
   ierr  = dbputqm(dbfile,trim(meshname), len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,db_f77null, & 
                   dims,ndims,db_float,db_noncollinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write non-collinear 2D quad mesh ', trim(meshname), &
                         ' to Silo file. Program execution stopped!'
      stop
   end if

   return

end subroutine mkQmesh2dNonCollinearStrain

! =======================================================================================================

subroutine mkQmesh2dNonCollinear(dbfile,meshname,outfileid,ndiv1,ndiv2,xcoord,ycoord)
   ! ===========================================================
   ! This routine creates a 2D quad mesh in an existing Silo
   ! file. The mesh is not collinear with the coordinate axes.
   !
   !    Input:
   !          dbfile            ID of the created Silo file
   !                            (integer value)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2       Number of cell divisions
   !                            along the two spatial
   !                            dimensions
   !          xcoord            x-coordinates of the nodes
   !          ycoord            y-coordinates of the nodes
   !
   !    Output:
   !          
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2
   !double precision, intent(in) :: xcoord((ndiv1+1)*(ndiv2+1)), ycoord((ndiv1+1)*(ndiv2+1))
   double precision, intent(in) :: xcoord(:), ycoord(:)
   character(*), intent(in)     :: meshname

   !double precision             :: xc(ndiv1+1,ndiv2+1), yc(ndiv1+1,ndiv2+1)
   real                         :: xc(ndiv1+1,ndiv2+1), yc(ndiv1+1,ndiv2+1)
   integer                      :: ierr, err, ndims, dims(2), dim1, dim2, i
 
   ! Write quad mesh data
   ! ====================
   ndims = 2
   dims  = (/ ndiv1+1, ndiv2+1 /)
   i     = 1
   do dim2=1,dims(2)
      do dim1=1,dims(1)
         xc(dim1,dim2) = xcoord(i)
         yc(dim1,dim2) = ycoord(i)
         i = i + 1
      end do
   end do
   
   ierr  = dbputqm(dbfile,trim(meshname),len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,db_f77null, &
                   dims,ndims,db_float,db_noncollinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write non-collinear 2D quad mesh ', trim(meshname), &
                         ' to Silo file. Program execution stopped!'
      stop
   end if

   return

end subroutine mkQmesh2dNonCollinear

! =======================================================================================================

subroutine mkQmesh3dNonCollinear(dbfile,meshname,outfileid,ndiv1,ndiv2,ndiv3,xcoord,ycoord,zcoord)
   ! ===========================================================
   ! This routine creates a 3D quad mesh in an existing Silo
   ! file. The mesh is not collinear with the coordinate axes.
   !
   !    Input:
   !          dbfile            ID of the created Silo file
   !                            (integer value)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2,ndiv3 Number of cell divisions
   !                            along the two spatial
   !                            dimensions
   !          xcoord            x-coordinates of the nodes
   !          ycoord            y-coordinates of the nodes
   !          zcoord            z-coordinates of the nodes
   !
   !    Output:
   !       
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3
   character(*), intent(in)     :: meshname
   double precision, intent(in) :: xcoord(:), ycoord(:), zcoord(:)
   !double precision, intent(in) :: xcoord((ndiv1+1)*(ndiv2+1)*(ndiv3+1))
   !double precision, intent(in) :: ycoord((ndiv1+1)*(ndiv2+1)*(ndiv3+1))
   !double precision, intent(in) :: zcoord((ndiv1+1)*(ndiv2+1)*(ndiv3+1))
   real                         :: xc(ndiv1+1,ndiv2+1,ndiv3+1)
   real                         :: yc(ndiv1+1,ndiv2+1,ndiv3+1)
   real                         :: zc(ndiv1+1,ndiv2+1,ndiv3+1)
   integer                      :: ierr, err, ndims, dims(3), dim1, dim2, dim3, i

   ! Write quad mesh data
   ! ====================
   ndims = 3
   dims  = (/ ndiv1+1, ndiv2+1, ndiv3+1 /)

   i = 1  
   do dim3=1,dims(3)
      do dim2=1,dims(2)
         do dim1=1,dims(1)
            xc(dim1,dim2,dim3) = xcoord(i)
            yc(dim1,dim2,dim3) = ycoord(i)
            zc(dim1,dim2,dim3) = zcoord(i)
            i = i + 1
         end do
      end do
   end do
   
   ierr  = dbputqm(dbfile,trim(meshname),len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,zc,dims, &
                   ndims,db_double,db_noncollinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write non-collinear 3D quad mesh ', trim(meshname), &
                         ' to Silo file. Program execution stopped!'
      stop
   end if

   return

end subroutine mkQmesh3dNonCollinear

! =======================================================================================================

subroutine mkQmesh3dNonCollinearStrain(dbfile,meshname,outfileid,ndiv1,ndiv2,ndiv3,xcoord,ycoord,zcoord,epsx,epsy,epsz)
   ! ===========================================================
   ! This routine creates a 3D quad mesh in an existing Silo
   ! file. The mesh is not collinear with the coordinate axes.
   !
   !    Input:
   !          dbfile            ID of the created Silo file
   !                            (integer value)
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          ndiv1,ndiv2,ndiv3 Number of cell divisions
   !                            along the two spatial
   !                            dimensions
   !          xcoord            x-coordinates of the nodes
   !          ycoord            y-coordinates of the nodes
   !          zcoord            z-coordinates of the nodes
   !          epsx              Strain component in the x-dir.
   !          epsy              Strain component in the y-dir.
   !          epsz              Strain component in the z-dir.
   !
   !    Output:
   !       
   !
   ! Last rev.:
   !    2010-09-01: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3
   double precision, intent(in) :: xcoord((ndiv1+1)*(ndiv2+1)*(ndiv3+1)), ycoord((ndiv1+1)*(ndiv2+1)*(ndiv3+1))
   double precision, intent(in) :: zcoord((ndiv1+1)*(ndiv2+1)*(ndiv3+1)), epsx, epsy, epsz
   real                         :: xc(ndiv1+1,ndiv2+1,ndiv3+1), yc(ndiv1+1,ndiv2+1,ndiv3+1), zc(ndiv1+1,ndiv2+1,ndiv3+1)
   character(*), intent(in)     :: meshname
   integer                      :: ierr, err, ndims, dims(3), dim1, dim2, dim3, i

   ! Write quad mesh data
   ! ====================
   ndims = 3
   dims  = (/ ndiv1+1, ndiv2+1, ndiv3+1 /)

   i = 1  
   do dim3=1,dims(3)
      do dim2=1,dims(2)
         do dim1=1,dims(1)
            xc(dim1,dim2,dim3) = xcoord(i)*epsx
            yc(dim1,dim2,dim3) = ycoord(i)*epsy
            zc(dim1,dim2,dim3) = zcoord(i)*epsz
            i = i + 1
         end do
      end do
   end do
   
   ierr  = dbputqm(dbfile,trim(meshname),len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,zc,dims,ndims, &
                   db_float,db_noncollinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write non-collinear 3D quad mesh ', trim(meshname), &
                         ' to Silo file. Program execution stopped!'
      stop
   end if


   return

end subroutine mkQmesh3dNonCollinearStrain

! =======================================================================================================

subroutine mkQmesh3dShiftedCollinearStrain(dbfile,meshname,outfileid,dim1,dim2,dim3,ndiv1,ndiv2,ndiv3, &
                                        start1,start2,start3,eps1,eps2,eps3)
   ! ===========================================================
   ! This routine creates a 3D quad mesh collinear with the
   ! coordinate axes. The mesh origin is set at the Cartesian
   ! coordinates (start1,start2,start3).
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          dim1,dim2,dim3    3D dimensions of the domain
   !          ndiv1,ndiv2,ndiv3 Number of element/cell
   !                            divisions along the three
   !                            spatial dimensions
   !          start1, start2,   The Cartesian coordinates
   !          start3            of the mesh origin
   !          eps1, eps2, eps3  Deformation ("strain") along
   !                            each coordinate axis
   !
   !    Output:
   !
   ! Last rev.:
   !    2011-04-12: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3
   double precision, intent(in) :: dim1, dim2, dim3, start1, start2, start3, eps1, eps2, eps3
   character(*), intent(in)     :: meshname
   integer                      :: ierr, err, ndims, dims(3), i
   double precision             :: xc(ndiv1+1), yc(ndiv2+1), zc(ndiv3+1), dx, dy, dz

   ! Write quad mesh data
   ! ====================
   ndims = 3
   dims  = (/ ndiv1+1, ndiv2+1, ndiv3+1 /)
   xc = start1
   yc = start2
   zc = start3
   dx = dim1/ndiv1
   dy = dim2/ndiv2
   dz = dim3/ndiv3
   do i=1,ndiv1
      xc(i+1) = (xc(i) + dx) 
   end do
   do i=1,ndiv2
      yc(i+1) = (yc(i) + dy) 
   end do
   do i=1,ndiv3
      zc(i+1) = (zc(i) + dz) 
   end do
   ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,xc*eps1,yc*eps2,zc*eps3,dims,ndims, &
                  db_double,db_collinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write shifted, deformed, collinear 3D quad mesh ', trim(meshname), &
                         ' to Silo file. Program execution stopped!'
      stop
   end if
   
   return

end subroutine mkQmesh3dShiftedCollinearStrain

! =======================================================================================================

subroutine mk3dQuadMultiMeshCollinear(dbfile,path,meshname,outfileid,nmesh)
   ! ===========================================================
   ! This routine creates a 3D multimesh in an existing Silo db-
   ! file. The mesh is collinear with the coordinate axes.
   !
   !    Input:
   !          dbfile            ID of the created Silo file
   !                            (integer value)
   !          path              Search path to the multimesh
   !                            files (character string)
   !          meshname          Name of the master mesh
   !                            (character string)
   !          outfileid         File-ID for program output
   !                            (integer value)
   !          nmesh             Number of meshes (integer value)
   !
   !    Output:
   !       
   !
   ! Last rev.:
   !    2010-12-20: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)      :: dbfile, nmesh, outfileid
   character(*), intent(in) :: meshname, path
   character(50)            :: meshnames(nmesh), buffer, ci
   integer                  :: ierr, err, lmeshnames(nmesh), meshtypes(nmesh), i, oldlen
   character, parameter     :: delimiter='/'

   ! Set length of mesh name strings
   oldlen = dbget2dstrlen()
   err    = dbset2dstrlen(50)

   ! Set meshnames and meshtypes
   do i=1,nmesh
      write(ci,'(I10)') i
      write(buffer,'(A10,A10)') 'subdomain.',adjustl(ci)
      meshnames(i) = trim(path) // delimiter // trim(buffer) // ':' // trim(meshname)
      !meshnames(i) = trim(buffer) // ':' // '/' // trim(meshname)
      lmeshnames(i) = len_trim(meshnames(i))
      meshtypes(i) = db_quad_rect

   end do

   ! Write multimesh objects
   err = dbputmmesh(dbfile, meshname, len_trim(meshname), nmesh, meshnames, lmeshnames, meshtypes, db_f77null, ierr)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not set multimesh. Program execution stopped!'
      stop
   end if
   ! Restore previous value for maximum string length
   err = dbset2dstrlen(oldlen)

   return

end subroutine mk3dQuadMultiMeshCollinear

! =======================================================================================================

subroutine mkQmesh3dShiftedCollinear(dbfile,meshname,outfileid,dim1,dim2,dim3,ndiv1,ndiv2,ndiv3,start1,start2,start3)
   ! ===========================================================
   ! This routine creates a 3D quad mesh collinear with the
   ! coordinate axes. The mesh origin is set at the Cartesian
   ! coordinates (start1,start2,start3).
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          dim1,dim2,dim3    3D dimensions of the domain
   !          ndiv1,ndiv2,ndiv3 Number of element/cell
   !                            divisions along the three
   !                            spatial dimensions
   !          start1, start2,   The Cartesian coordinates
   !          start3            of the mesh origin
   !
   !    Output:
   !
   ! Last rev.:
   !    2010-12-15: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, outfileid, ndiv1, ndiv2, ndiv3
   double precision, intent(in) :: dim1, dim2, dim3, start1, start2, start3
   character(*), intent(in)     :: meshname
   integer                      :: ierr, err, ndims, dims(3), i
   double precision             :: xc(ndiv1+1), yc(ndiv2+1), zc(ndiv3+1), dx, dy, dz

   ! Write quad mesh data
   ! ====================
   ndims = 3
   dims  = (/ ndiv1+1, ndiv2+1, ndiv3+1 /)
   xc = start1
   yc = start2
   zc = start3
   dx = dim1/ndiv1
   dy = dim2/ndiv2
   dz = dim3/ndiv3
   do i=1,ndiv1
      xc(i+1) = xc(i) + dx
   end do
   do i=1,ndiv2
      yc(i+1) = yc(i) + dy
   end do
   do i=1,ndiv3
      zc(i+1) = zc(i) + dz
   end do
   ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,zc,dims,ndims, &
                  db_double,db_collinear,db_f77null,err)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not write shifted collinear 3D quad mesh ', trim(meshname), &
                         ' to Silo file. Program execution stopped!'
      stop
   end if
   
   return

end subroutine mkQmesh3dShiftedCollinear

! =======================================================================================================

subroutine setMultiVarQuad(dbfile,varname,path,outfileid,nvars)
   ! ===========================================================
   ! This routine writes integer multivars to an existing
   ! master Silo db-file with a quad mesh.
   !
   !    Input:
   !          dbfile            ID of the existing Silo file
   !                            (integer value)
   !          path              Search path to the multimesh
   !                            files (character string)
   !          meshname          Name of the master mesh
   !                            (character string)
   !          outfileid         File-ID for program output
   !                            (integer value)
   !          nmesh             Number of meshes (integer value)
   !
   !    Output:
   !       
   !
   ! Last rev.:
   !    2010-12-20: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)      :: dbfile, nvars, outfileid
   character(*), intent(in) :: varname, path
   character(50)            :: varnames(nvars), buffer, ci
   integer                  :: ierr, err, lvarnames(nvars), vartypes(nvars), i, oldlen
   character, parameter     :: delimiter='/'

   ! Set length of varname strings
   oldlen = dbget2dstrlen()
   err    = dbset2dstrlen(50)

   ! Set varnames and vartypes
   do i=1,nvars
      write(ci,'(I10)') i
      write(buffer,'(A10,A10)') 'subdomain.',adjustl(ci)
      varnames(i) = trim(path) // delimiter // trim(buffer) // ':' // trim(varname)
      lvarnames(i) = len_trim(varnames(i))
      vartypes(i) = db_quadvar

   end do

   ! Write multimesh objects
   err = dbputmvar(dbfile, varname, len_trim(varname), nvars, varnames, lvarnames, vartypes, db_f77null, ierr)
   if (ierr.eq.-1) then
      write(outfileid,*) 'Could not set multivar. Program execution stopped!'
      stop
   end if
   ! Restore previous value for maximum string length
   err = dbset2dstrlen(oldlen)

   return

end subroutine setMultiVarQuad

! =======================================================================================================

subroutine mkQmesh(dbfile,meshname,outfileid,coord,coll)
   ! ===========================================================
   ! This routine creates a 3D quad mesh collinear with the
   ! coordinate axes. The mesh origin is set at the Cartesian
   ! coordinates (start1,start2,start3).
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          outfileid         File-ID for program output 
   !          coord             Node coordinates [3 x nnodes]
   !             [x1 x2 ... xn
   !              y1 y2 ... yn
   !              z1 z2 ... zn]
   !          coll = 0          Collinear mesh
   !                 1          Non-collinear mesh
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-05: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)           :: dbfile, outfileid, coll
   double precision, intent(in)  :: coord(:,:)
   character(*), intent(in)      :: meshname
   integer                       :: ierr, err, ndims, dim, dims(size(coord,1)), i, dim1, dim2, dim3
   double precision, allocatable :: xc(:,:,:), yc(:,:,:), zc(:,:,:), xc2(:,:), yc2(:,:)

   ! Write quad mesh data
   ! ====================
   ndims = size(coord,1)
   dim   = size(coord,2)

   if (ndims.eq.2) then ! 2D mesh
   ! -----------------------------------------------------------------------------------------------------     
     
      dims  = (/ dim, dim, dim /)
   
      select case (coll)
   
         case(0) ! Collinear mesh
         
            ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,coord(1,:),coord(:,2),coord(:,3),dims,ndims, &
                           db_double,db_collinear,db_f77null,err)
                          
         case(1) ! Non-collinear mesh
        
            allocate(xc(dim,dim,dim),yc(dim,dim,dim),zc(dim,dim,dim))
            i = 1  
            do dim3=1,dim
               do dim2=1,dim
                  do dim1=1,dim
                     xc(dim1,dim2,dim3) = coord(1,i)
                     yc(dim1,dim2,dim3) = coord(2,i)
                     zc(dim1,dim2,dim3) = coord(3,i)
                    i = i + 1
                  end do
               end do
            end do
         
            ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,xc,yc,zc,dims,ndims, &
                           db_double,db_noncollinear,db_f77null,err)
         
            deallocate(xc, yc, zc)
         
         case default
         
            write(outfileid,*) 'Wrong value of variable "coll" (0 or 1) in subroutine "mkQmesh". Program execution stopped!'
            stop
         
      end select
                    
      if (ierr.eq.-1) then
         write(outfileid,*) 'Could not write collinear quad mesh ', trim(meshname), &
                            ' to Silo file. Error in subroutine "mkQmesh". Program execution stopped!'
         stop
      end if
      
   else if (ndims.eq.3) then ! 3D mesh
   ! -----------------------------------------------------------------------------------------------------
     
      dims  = (/ dim, dim /)
     
      select case (coll)
   
         case(0) ! Collinear mesh

            ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,coord(1,:),coord(2,:),db_f77null, &                
                           dims,ndims,db_double,db_collinear,db_f77null,err)

         case(1) ! Non-collinear mesh
        
            allocate(xc2(dim,dim),yc2(dim,dim))
            i = 1  
            do dim2=1,dim
               do dim1=1,dim
                  xc2(dim1,dim2) = coord(1,i)
                  yc2(dim1,dim2) = coord(2,i)
                  i = i + 1
               end do
            end do
         
            ierr = dbputqm(dbfile,meshname,len_trim(meshname),'xc',2,'yc',2,'zc',2,xc2,yc2,db_f77null,dims,ndims, &
                           db_double,db_noncollinear,db_f77null,err)
         
            deallocate(xc, yc)
         
         case default
        
            write(outfileid,*) 'Wrong value of variable "coll" (0 or 1) in subroutine "mkQmesh2d". Program execution stopped!'
            stop
         
      end select
                    
      if (ierr.eq.-1) then
         write(outfileid,*) 'Could not write collinear 2D quad mesh ', trim(meshname), &
                            ' to Silo file. Error in subroutine "mkQmesh2d". Program execution stopped!'
         stop
      end if

     
   else
   ! -----------------------------------------------------------------------------------------------------     
   
         write(outfileid,*) 'Wrong dimension of "coord" in subroutine "mkQmesh". Program execution stopped!'
         stop

   end if
       
   return

end subroutine mkQmesh

! =======================================================================================================

subroutine plotSolidMesh(dbfile,meshname,coord,enode)
   ! ===========================================================
   ! This routine creates an unstructured 2D or 3D mesh. The
   ! following element types are supported:
   !    2D: 3-node triangles and 4-node quads
   !    3D: 4-node tetrahedrons and 8-node bricks
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          coord             Node coordinates [2 or 3 x nnodes]
   !             [x1 x2 ... xn
   !              y1 y2 ... yn
   !             (z1 z2 ... zn)]
   !          enode             Element nodes [nelnodes x nelem]
   !             [n1 n1 ...
   !              n2 n2 ...
   !              n3 n3 ...
   !              :  :
   !              n8 n8 ...]
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile, enode(:,:)
   double precision, intent(in) :: coord(:,:)
   character(*), intent(in)     :: meshname
   integer                      :: ierr, err, ndims, nodelist(size(enode,1)*size(enode,2)), nenodes, i, j
   integer                      :: lnodelist, shapesize, shapecounts, nshapetypes, nzones, nnodes, shapetype
     
   character(92)                :: errmsg_1 = 'Could not write connectivity information in subroutine "mkUmesh". Program execution stopped!'
   character(85)                :: errmsg_2 = 'Could not write unstructured mesh in subroutine "mkUmesh". Program execution stopped!'
   character(81)                :: errmsg_3 = 'Wrong number of element nodes in subroutine "mkUmesh". Program execution stopped!'
   character(96)                :: errmsg_4 = 'Wrong number of rows in the coordinate array in subroutine "mkUmesh". Program execution stopped!'
   
   nzones      = size(enode,2)    ! Number of zones or elements
   shapecounts = nzones           ! Number of zone shapes
   lnodelist   = size(nodelist,1) ! Element (zone) nodes
   nshapetypes = 1                ! Only one type of elements present
   nnodes      = size(coord,2)    ! Number of nodes
   
   select case (size(coord,1))    ! Check if the current mesh is 2D or 3D

! ===================================================================================================================   
      case (2) ! 2D mesh
! ===================================================================================================================   
            
         ndims = 2                ! 2-dimensional mesh
            
         select case (size(enode,1))
   
            case (3) ! Triangular elements
        
               shapetype = db_zonetype_triangle
               nenodes   = 3
               shapesize = 3
         
            case (4) ! Quad elements
        
               shapetype = db_zonetype_quad
               nenodes   = 4
               shapesize = 4
               
            case default ! Wrong number of nodes
        
               write(*,*) errmsg_3
               stop
         
         end select

         ! Create Silo-formatted node list
         ! -------------------------------
         j = 1
         do i=1,nzones
            nodelist(j:nenodes*i) = enode(:,i)
            j = nenodes*i + 1
         end do

         ! Write connectivity information
         ! ------------------------------
         ierr = dbputzl2(dbfile,'zonelist',8,nzones,ndims,nodelist,lnodelist,1,0,0, &
                         shapetype,shapesize,shapecounts,nshapetypes,db_f77null,err)
                    
         if (ierr.eq.-1) then
            write(*,*) errmsg_1
            stop
         end if
        
         ! Write unstructured mesh
         ! -----------------------
         ierr = dbputum(dbfile,trim(meshname),len_trim(meshname),ndims,coord(1,:),coord(2,:),db_f77null,'xc',2,'yc',2,db_f77null,0, &
                        db_double,nnodes,nzones,'zonelist',8,db_f77null,0,db_f77null,err)

         if (ierr.eq.-1) then
            write(*,*) errmsg_2
            stop
         end if

! ===================================================================================================================         
      case (3) ! 3D mesh
! ===================================================================================================================

         ndims = 3                ! 3-dimensional mesh
        
         select case (size(enode,1))
   
            case (4) ! Tetrahedral elements
        
               shapetype = db_zonetype_tet
               nenodes   = 4
               shapesize = 4
         
            case (8) ! Brick elements
        
               shapetype = db_zonetype_hex
               nenodes   = 8
               shapesize = 8
               
            case default ! Wrong number of nodes
        
               write(*,*) errmsg_3
               stop
         
         end select
        
         ! Create Silo-formatted node list
         ! -------------------------------
         j = 1
         do i=1,nzones
            nodelist(j:nenodes*i) = enode(:,i)
            j = nenodes*i + 1
         end do

         ! Write connectivity information
         ! ------------------------------
         ierr = dbputzl2(dbfile,'zonelist',8,nzones,ndims,nodelist,lnodelist,1,0,0, &
                         shapetype,shapesize,shapecounts,nshapetypes,db_f77null,err)
                    
         if (ierr.eq.-1) then
            write(*,*) errmsg_1
            stop
         end if
        
         ! Write unstructured mesh
         ! -----------------------
         ierr = dbputum(dbfile,trim(meshname),len_trim(meshname),ndims,coord(1,:),coord(2,:),coord(3,:),'xc',2,'yc',2,'zc',2, &
                        db_double,nnodes,nzones,'zonelist',8,db_f77null,0,db_f77null,err)

         if (ierr.eq.-1) then
            write(*,*) errmsg_2
            stop
         end if
        
      case default
        
         write(*,*) errmsg_4
         stop
         
   end select
   
   return

end subroutine plotSolidMesh

! =======================================================================================================

subroutine setUmeshDVar(dbfile,meshname,var,varname,vartype)
   ! ===========================================================
   ! This routine writes scalar values to an unstructured 2D or
   ! 3D mesh. Nodal or zonal double precision values are
   ! supported.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          var               Node or element variable
   !                            [nnodes x 1] or [nel x 1]
   !          varname           Variable name
   !          vartype           Type of variable:
   !                            'zonal' or 'nodal'
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile
   double precision, intent(in) :: var(:)
   character(*), intent(in)     :: meshname, varname, vartype
   integer                      :: ierr, err
   
   if (vartype.eq.'zonal') then
     
      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_double, db_zonecent, db_f77null, ierr)
      
   else if (vartype.eq.'nodal') then

      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_double, db_nodecent, db_f77null, ierr)
     
   else
     
      write(*,*) 'Wrong variable "vartype" in subroutine setUmeshDVar". Program execution stopped!'
      stop
      
   end if
   
   return

end subroutine setUmeshDVar

! =======================================================================================================

subroutine setUmeshIVar(dbfile,meshname,var,varname,vartype)
   ! ===========================================================
   ! This routine writes scalar values to an unstructured 2D or
   ! 3D mesh. Nodal or zonal integer values are
   ! supported.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          var               Node or element variable
   !                            [nnodes x 1] or [nel x 1]
   !          varname           Variable name
   !          vartype           Type of variable:
   !                            'zonal' or 'nodal'
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-09: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)      :: dbfile
   integer, intent(in)      :: var(:)
   character(*), intent(in) :: meshname, varname, vartype
   integer                  :: ierr, err
   
   if (vartype.eq.'zonal') then
     
      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_int, db_zonecent, db_f77null, ierr)
      
   else if (vartype.eq.'nodal') then

      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_int, db_nodecent, db_f77null, ierr)
     
   else
     
      write(*,*) 'Wrong variable "vartype" in subroutine setUmeshIVar". Program execution stopped!'
      stop
      
   end if
   
   return

end subroutine setUmeshIVar

! ===================================================================================================================   
            
subroutine setUmeshDVar2(dbfile,meshname,var,varname,vartype,varunit)
   ! ===========================================================
   ! This routine writes scalar values to an unstructured 2D or
   ! 3D mesh. Nodal or zonal double precision values are
   ! supported.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          var               Node or element variable
   !                            [nnodes x 1] or [nel x 1]
   !          varname           Variable name
   !          vartype           Type of variable:
   !                            'zonal' or 'nodal'
   !          varunit           Variable unit
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-11: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile
   double precision, intent(in) :: var(:)
   character(*), intent(in)     :: meshname, varname, vartype, varunit
   integer                      :: ierr, err
   integer                      :: optlistid
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(1,optlistid) 
   err = dbaddcopt(optlistid, dbopt_units, varunit, len_trim(varunit))
   
   if (vartype.eq.'zonal') then
     
      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_double, db_zonecent, optlistid, ierr)
      
   else if (vartype.eq.'nodal') then

      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_double, db_nodecent, optlistid, ierr)
     
   else
     
      write(*,*) 'Wrong variable "vartype" in subroutine setUmeshDVar". Program execution stopped!'
      stop
      
   end if
   
   ! Free options list
   ! =================
   err = dbfreeoptlist(optlistid)
   
   return

end subroutine setUmeshDVar2

! =======================================================================================================

subroutine setUmeshIVar2(dbfile,meshname,var,varname,vartype,varunit)
   ! ===========================================================
   ! This routine writes scalar values to an unstructured 2D or
   ! 3D mesh. Nodal or zonal integer values are
   ! supported.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          var               Node or element variable
   !                            [nnodes x 1] or [nel x 1]
   !          varname           Variable name
   !          vartype           Type of variable:
   !                            'zonal' or 'nodal'
   !          varunit           Variable unit
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-11: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)      :: dbfile
   integer, intent(in)      :: var(:)
   character(*), intent(in) :: meshname, varname, vartype, varunit
   integer                  :: ierr, err
   integer                  :: optlistid
   
   ! Create option list
   ! ==================
   err = dbmkoptlist(1,optlistid) 
   err = dbaddcopt(optlistid, dbopt_units, varunit, len_trim(varunit))
   
   if (vartype.eq.'zonal') then
     
      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_int, db_zonecent, optlistid, ierr)
      
   else if (vartype.eq.'nodal') then

      err = dbputuv1(dbfile, trim(varname), len_trim(varname), trim(meshname), len_trim(meshname), var, &
                     size(var,1), db_f77null, 0, db_int, db_nodecent, optlistid, ierr)
     
   else
     
      write(*,*) 'Wrong variable "vartype" in subroutine setUmeshIVar". Program execution stopped!'
      stop
      
   end if
   
   ! Free options list
   ! =================
   err = dbfreeoptlist(optlistid)
   
   return

end subroutine setUmeshIVar2

! ===================================================================================================================   

subroutine plotPointMesh(dbfile,meshname,coord)
   ! ===========================================================
   ! This routine writes a 2D or 3D point mesh to an existing
   ! Silo file.
   !
   !    Input:
   !          dbfile            Pointer to a Silo file
   !          meshname          Name of the mesh
   !          coord             Point coordinates [2 or 3 x nnodes]
   !             [x1 x2 ... xn
   !              y1 y2 ... yn
   !              (z1 z2 ... zn)]
   !
   !    Output:
   !          Data is written to an existing Silo file
   !
   ! Last rev.:
   !    2011-05-11: H. Hallberg
   ! ===========================================================

   implicit none

   integer, intent(in)          :: dbfile
   character(*), intent(in)     :: meshname
   double precision, intent(in) :: coord(:,:)
   integer                      :: ierr, err
   integer                      :: npoints, ndims
   
   ndims   = size(coord,1)
   npoints = size(coord,2)

   select case (ndims)
   
      case (2)
        
         ierr = dbputpm(dbfile,trim(meshname),len_trim(meshname),ndims,coord(1,:),coord(2,:),db_f77null,npoints, &
                        db_double,db_f77null,err)
         if (ierr.eq.-1) then
            write(*,*) 'Could not write 2D point mesh in subroutine "plotPoints". Program execution stopped!'
            stop
         end if
         
      case(3)
        
         ierr = dbputpm(dbfile,trim(meshname),len_trim(meshname),ndims,coord(1,:),coord(2,:),coord(3,:),npoints, &
                        db_double,db_f77null,err)
         if (ierr.eq.-1) then
            write(*,*) 'Could not write 3D point mesh in subroutine "plotPoints". Program execution stopped!'
            stop
         end if
      
      case default
        
         write(*,*) 'Wrong number of rows in "coord" in subroutine "plotPoints". Program execution stopped!'
         stop
         
   end select
   
   return

end subroutine plotPointMesh

! ===================================================================================================================   

end module wrt2silo
