module ioGMSH

  use structure

  implicit none

  integer, dimension(:), allocatable :: bocoToFlag_, perioToFlag_
  public readMesh, closeMesh

contains

  subroutine readMeshGMSH

    implicit none

    ! loc
    integer :: ios
    integer :: n_per_joins

    ! Temporay storage for periodic joins
    type (strJoinBoundaryCondition), dimension(:), allocatable :: per_joins ! Intern joins


    open(global%indexfile, file = global%gridFile, & 
         status = "old", action = "read", iostat = ios)

    if (ios .ne. 0) then
      print *,"File "//trim(adjustl(global%gridFile))//" not found!"
      stop
    endif

    print *,"_____________________________________________________________"
    print *,"                                                             "
    print *,"                             Mesh                            "
    print *,"_____________________________________________________________"
    print *,"==> Read and treat the mesh                                  "

    call readHeader

    call readFamily(n_per_joins, per_joins)

    call readDimension

    call readElements

  end subroutine readMeshGMSH



  !=============================================================================!
  ! ReadHeader
  !=============================================================================!
  subroutine readHeader

    implicit none

    ! local variables
    integer                                :: dummy1, dummy2
    real                                   :: version

    ! Reading the header of the mesh
    ! In the header, we are only interested in the version number
    read(global%indexfile,*)
    read(global%indexfile,*) version, dummy1, dummy2

    if (version .ne. 2.2) then
      call terminate("IoGMSH", "Bad version of GMSH: use the format 2.2")
    endif

    read(global%indexfile,*)

  end subroutine readHeader

  !=============================================================================!
  ! readFamily
  !=============================================================================!

  subroutine readFamily(n_per_joins, per_joins)

    implicit none


    integer :: ios
    integer :: nFams, dim, flag
    character (len = 100) :: gmshFamilyName
    type(strGMSH_Family), dimension(:), allocatable :: GMSH_Family
        integer, dimension(:), allocatable :: tmpBocoIdGmsh, tmpPerioIdGmsh, tmpBocoType


    call gotoSection(global%indexfile, '$PhysicalNames')
    
    read(global%indexfile, *, iostat=ios) nFams
    if (ios /= 0) then
       rewind(global%indexfile)
       return
    end if

    allocate(GMSH_Family(nFams - 1))
    allocate(bocoToFlag_(InputMarkers%nMarkers + InputPerMarkers%nPerMarkers))
    allocate(perioToFlag_(InputPerMarkers%nPerMarkers))))
    allocate(tmpBocoIdGmsh(InputMarkers%nMarkers + InputPerMarkers%nPerMarkers))
    allocate(tmpBocoType(InputMarkers%nMarkers + InputPerMarkers%nPerMarkers))
    allocate(tmpPerioIdGmsh(InputPerMarkers%nPerMarkers))


    nBocos = 0
    n_per_joins = 0

    do iFam_Gmsh = 1, nFams
      read(global%indexfile, *) dim, flag, gmshFamilyName

      ! Find corresponding BOCO
      do iMarker = 1, InputMarkers%nMarkers
        if (trim(adjustl(InputMarkers%name(iMarker))) .eq. trim(adjustl((gmshFamilyName)))) then
          nBocos = nBocos + 1
          !  tmpBocoType(nBocos) = InputMarkers%type(iMarker)
          !  tmpBocoIdGmsh(nBocos) = iFam_Gmsh

          GMSH_Family(iFam_Gmsh)%name     = trim(adjustl((gmshFamilyName)))
          GMSH_Family(iFam_Gmsh)%flag     = flag
          GMSH_Family(iFam_Gmsh)%idMarker = iMarker
          bocoToFlag_(nBocos) = flag
          exit
        end if
      end do

      ! Or periodic condition
      do iMarker = 1, InputPerMarkers%nPerMarkers
        if (trim(adjustl(InputPerMarkers%joinName(iMarker))) .eq. trim(adjustl(gmshFamilyName))) then
          nBocos = nBocos + 1
          tmpBocoType(nBocos) = PERIODIC
          tmpBocoIdGmsh(nBocos) = iFam_Gmsh

          bocoToFlag_(nBocos) = flag

          GMSH_Family(iFam_Gmsh)%name     = trim_ws(gmshFamilyName)
          GMSH_Family(iFam_Gmsh)%flag     = flag
          GMSH_Family(iFam_Gmsh)%idMarker = iMarker

          ! NB: we want join order to be the same as periodic marker order
          n_per_joins = n_per_joins + 1
          perioToFlag_(iMarker) = flag
          tmpPerioIdGmsh(iMarker) = iFam_Gmsh
          exit
        end if
      end do

    end do

    ! Allocate bocos
    allocate(Boundary(nBocos))

    Boundary(1:nBocos)%bocotype = tmpBocoType(1:nBocos)
    Boundary(1:nBocos)%idGmsh = tmpBocoIdGmsh(1:nBocos)

    ! Allocate periodic
    allocate(per_joins(n_per_joins))

    per_joins(1:n_per_joins)%idGmsh = tmpPerioIdGmsh(1:n_per_joins)

    rewind(global%indexfile)

  end subroutine readFamily


  !=============================================================================!
  ! gotoSection
  !=============================================================================!

  subroutine gotoSection(file_idx, section_name)
    integer, intent(in)               :: file_idx
    character (len = *), intent(in)   :: section_name
    integer               :: ios, str_length
    character (len = 100) :: line

    str_length = len(trim_ws(section_name))
    line = ' '
    do while(adjustl(line(1:str_length)) /= section_name)
       read (file_idx, *, iostat=ios) line
       if(ios /= 0) then
          !! if(zone%rank .eq. 0) print *, section_name//" not found"
          if ( section_name == "$PhysicalNames" ) then
            rewind(file_idx)
            return           
          else
            call terminate( "gotoSection", "In IoGMSH.f90, "//section_name//" not found" )
          endif
       end if
    end do

  end subroutine gotoSection



  !=============================================================================!
  ! readDimension
  !=============================================================================!

  subroutine readDimension

    implicit none

    integer :: ierr, nElt, nHexa, nQuad, eID, iType, k, ios, is3D=0

    ierr = 0
    nQuad = 0
    nHexa = 0
    Mesh%max_vpc = 0

    call gotoSection(global%indexfile, '$Nodes')
    read(global%indexfile,*) Mesh%n_vert_tot

    call gotoSection(global%indexfile, '$Elements')
    read(global%indexfile,*) nElt

    do k = 1, nElt
      read (global%indexfile, *)  eID, iType
      !Only hexa and quads for now
      if (iType .eq. GMSH_HEXA) then
        nHexa = nHexa + 1
        Mesh%max_vpc = max(Mesh%max_vpc, 8)
      else if (itype .eq. GMSH_HEXA_27N) then
        nHexa = nHexa + 1
        Mesh%max_vpc = max(Mesh%max_vpc, 27)
      else if (iType .eq. GMSH_QUAD) then
        nQuad = nQuad + 1
        Mesh%max_vpc = max(Mesh%max_vpc, 4)
      else if (itype .eq. GMSH_QUAD_9N) then
        nQuad = nQuad + 1
        Mesh%max_vpc = max(Mesh%max_vpc, 9)
      end if
    end do

    if(nHexa .eq. 0) then
       is3D = 0
       Mesh%n_cell_tot = nQuad
    else
       is3D = 1
       Mesh%n_cell_tot = nHexa
    end if

    rewind(unit = global%indexFile, iostat = ios)

    if(is3D .eq. 1) then
      print *,"The mesh contains 3D elements. Not allowed!"; stop
    end if

  end subroutine readDimension

  end module ioGMSH
