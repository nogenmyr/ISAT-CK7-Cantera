!---------------------------------------------------------------------!
! OWNER: Cornell University                                           !
! COPYRIGHT: © 2012, Cornell University                               !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-CU.txt' file included in the ISAT-CK7     !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

!
! STREAMS_MOD: MODULE TO INTERFACE THE NEW STREAMS FILE FORMAT FOR ISAT-CK
!
! GUIDE FOR USING STREAMS MODULE
! 
! Haifeng Wang
! Cornell University
! March, 2010
!
! INTRODUCTION
!
! The streams module interprets the new format of the input streams file 
! (streams.in), and provides the interface for ISAT-CK to get the streams data. 
! Two types of interfaces are provided for ISAT-CK to interact with this 
! module: public functions and public data.
!
! PUBLIC FUNCTIONS
!
! The following public functions are provided for ISAT-CK to call. 
!
! 1. subroutine streams_init(Unit,err,NAMES,NNAMES)
!
!     Initialize the streams module so that all the input from streams.in 
!     are interpreted and stored (in public data StrmDat). The optional 
!     arguments NAMES and NNAMES must be provided for all modes except 
!     modeci=1 (CONST_PROPERTY). 
!
! 2. function   streams_status() 
!
!     check the status of the module. It returns .true. if the module is 
!     initialized successfully after calling streams_init; otherwise it 
!     returns .false.. It is strongly recommended to always check the status 
!     before using any data or calling functions in this module.
!
! 3. subroutine streams_free()
!	
!     Free the data stored in this module after they are not needed anymore. 
!     After this function is called, streams_status returns .false.	
!
! 4. function   streams_index(name,err)
!
!	  Return the stream index given the stream name
!
! 5. function   streams_name(istream,err)
!
!     Return the stream name given the stream index
!
! 6. subroutine streams_ErrMsg(line,error)
!	
!     Retrieve the error message and the line in streams.in that contains the error.
!
! 7. function   streams_format(Unit,err)
!
!     Special function to return the format of streams.in without the 
!     initialization of the module, i.e. before calling streams_init
!
! 8. function   streams_modeci(Unit,err)
!
!     Special function to return the mode of ISAT-CK (MODECI) without the 
!     initialization of the module, i.e. before calling streams_init
!
! PUBLIC DATA
!
! The streams data are stored in the derived data type variable StrmDat. The 
! attribute of StrmDat is public, so it can be accessed by ISAT-CK directly.
! The status of the module should be checked before it is used (by calling 
! streams_status) to ensure the proper initialization of StrmDat. StrmDat contains
! the following components
!
! StrmDat%iformat          format of streams.in
!                          applicable to all modeci
!
! StrmDat%ModeCI           mode of ISAT-CK
!                          = 1 or CONST_PROPERTY     - inert, constant-density flow
!                          = 2 or MIXTURE_FRACTION   - mixture fraction formulation
!                          = 3 or PROGRESS_VARIABLE  - reaction progress variable formulation 
!                          = 4                       - not supported yet
!                          = 5                       - not supported yet
!                          = 6 or DIRECT_INTEGRATION - CHEMKIN/direct integration
!                          = 7 or ISAT_DI            - CHEMKIN/ISAT
!                          = 8 or RCCE               - CHEMKIN/ISAT/RCCE
!                          = 9 or ICE_PIC            - CHEMKIN/ISAT/ICE_PIC
!
! StrmDat%density          density value
!                          applicable to mode {1}
!
! StrmDat%timescale        reaction time scale
!                          applicable to mode {3}
!
! StrmDat%nStream          number of streams
!                          applicable to modes {2, 3, 6, 7, 8, and 9}
!
! StrmDat%nC               size of each stream
!                          applicable to modes {2, 3, 6, 7, 8, and 9}
!                          StrmDat%nC = NNAMES which is from streams_init
!
! StrmDat%Cvalue(:,:)      store the streams data 
!                          applicable to modes {2, 3, 6, 7, 8, and 9}
!                          one stream per row
!                          size=(nC,nStream)
!
! StrmDat%kstrm(:)         integer flag to indicate equilibrium state
!                          applicable to modes {6, 7, 8, and 9}
!                          = 1: the stream composition is that specified
!                          = 2: the stream composition is the equilibrium state
!                               (given the same P, enthalpy)
!                          = 3: the stream composition is the equilibrium state
!                               (given the same P, temperature)
!                          size=(nStream)
! 
! StrmDat%kmole(:)         integer flag to indicate the unit for the input species values
!                          applicable to modes {6, 7, 8, and 9}
!                          = 1: relative moles
!                          = 0: relative mass  
!                          size=(nStream)
!
! StrmDat%nSp              number of species (or represented species)
!                          applicable to modes {2, 3, 8, and 9}
!                          for modes={2,3}, it is the number of species in the flame table
!                          for modes={8,9}, it is the number of represented species
!
! StrmDat%Spnames(:)       names of species (or represented species)
!                          applicable to modes {2, 3, 8, and 9}
!                          for modes={2,3}, it is the names of species in the flame table
!                          for modes={8,9}, it is the names of represented species
!
! StrmDat%nTable           flame table size (number of table rows)
!                          applicable to modes {2, 3}
!
! StrmDat%Table(:,:)       flame table data
!                          applicable to modes 2, 3
!                          size=(nSp+4,nTable)
!                          each row contains: f, dens, p, T, species
! 
! MISCELLANEOUS
!
! 1. Specify NAMES,NNAMES for streams_init
! 
!     The optional arguments NAMES, NNAMES are required for all modes except 
!     (1, CONST_PROPERTY). NAMES contains all possible input variables that can
!     possibly present in the streams input (between STREAMS BEGIN/END in 
!     streams.in). For modes 2 or 3, NAMES should contain at least 'f'. 
!     For modes 6, 7, 8, and 9, NAMES should contain at least 'P', 'T', and one 
!     species name. The values for each input variable (e.g., 'P', 'T') in 
!     streams.in are stored in StrmDat%Cvalue in the corresponding location 
!     specified by NAMES. For those variables not having the corresponding input
!     from streams.in, their values in StrmDat%Cvalue are by default zero.
!
! 2. Error handling
!     
!     Most public functions in this module return the error status (err). After 
!     calling those functions, the error status should be checked in ISAT-CK
!     to ensure appropriate usage of this module. If error occurs (err=.true.),
!     ISAT-CK should abort the execution and display the error message in its
!     consistent style. The error message and the line in streams.in containing 
!     the error can be retrieved from streams_ErrMsg(line,error)
!     
! 3. Template for streams.in (mode=1)
!
!     !----------------------------------------------------------------------------!  
!     !                                                                            !
!     !      Template for streams.in file (new format)                             !
!     !      (for mode 1)                                                          !
!     !      HW, March, 2010                                                       !
!     !                                                                            !
!     !      Note: * any part beginning with "!" is a comment                      !
!     !            * empty lines are ignored                                       !
!     !            * keywords are not case-sensitive                               !
!     !                                                                            !
!     !----------------------------------------------------------------------------!
!     !                                                                   
!     ! specify MODECI: the mode of thermochemistry to be used in ISAT
!     !        = 1 or CONST_PROPERTY     - inert, constant-density flow
!     !        = 2 or MIXTURE_FRACTION   - mixture fraction formulation
!     !        = 3 or PROGRESS_VARIABLE  - reaction progress variable formulation 
!     !        = 4                       - not supported yet
!     !        = 5                       - not supported yet
!     !        = 6 or DIRECT_INTEGRATION - CHEMKIN/direct integration
!     !        = 7 or ISAT_DI            - CHEMKIN/ISAT
!     !        = 8 or RCCE               - CHEMKIN/ISAT/RCCE
!     !        = 9 or ICE_PIC            - CHEMKIN/ISAT/ICE_PIC
!     !              
!     MODECI      1
!     !
!     ! specify density
!     !
!     DENSITY     1.2
!     !- end of streams.in file
!     
! 4. Template for streams.in (modes=2,3)
!
!     !----------------------------------------------------------------------------!
!     !                                                                            !
!     !      Template for streams.in file (new format)                             !
!     !      (for modes 2 and 3)                                                   !
!     !      HW, March, 2010                                                       !
!     !                                                                            !
!     !      Note: * any part beginning with "!" is a comment                      !
!     !            * empty lines are ignored                                       !
!     !            * keywords are not case-sensitive                               !
!     !                                                                            !
!     !----------------------------------------------------------------------------!
!     !                                                                   
!     ! specify MODECI: the mode of thermochemistry to be used in ISAT
!     !        = 1 or CONST_PROPERTY     - inert, constant-density flow
!     !        = 2 or MIXTURE_FRACTION   - mixture fraction formulation
!     !        = 3 or PROGRESS_VARIABLE  - reaction progress variable formulation 
!     !        = 4                       - not supported yet
!     !        = 5                       - not supported yet
!     !        = 6 or DIRECT_INTEGRATION - CHEMKIN/direct integration
!     !        = 7 or ISAT_DI            - CHEMKIN/ISAT
!     !        = 8 or RCCE               - CHEMKIN/ISAT/RCCE
!     !        = 9 or ICE_PIC            - CHEMKIN/ISAT/ICE_PIC
!     !            
!     MODECI      3
!     !
!     ! specify timescale (for mode=3 or PROGRESS_VARIABLE only)
!     !
!     TIMESCALE  0.3
!     !
!     !      Specify streams.
!     !      - Each stream begins with 'STREAM BEGIN' and ends with 'STREAM END'.
!     !      - Each stream should have a unique name, e.g., "FUEL", "PILOT", "COFLOW"
!     !      - f  : mixture fraction or progress variables
!     !
!     ! first stream: PILOT 
!     STREAM BEGIN PILOT 
!         f          0.27               ! progress variable or mixture fraction
!     STREAM END PILOT
!     ! second stream: COFLOW
!     STREAM BEGIN COFLOW 
!         f          0.
!     STREAM END COFLOW
!     ! third stream: FUEL
!     STREAM BEGIN FUEL
!         f          1.
!     STREAM END FUEL
!     !
!     ! names of the species in the flame table
!     !
!     TABLE_VARIABLES BEGIN
!       H2 O2
!       H2O
!     TABLE_VARIABLES END
!     !
!     ! flame table data
!     ! note: no comment lines or empty lines are allowed when specifying the table data
!     ! f, dens, p, T, species
!     !
!     TABLE_DATA BEGIN
!      0.0 1.0 1.0 300.0 1.0 0.0 0.0
!      0.2 1.0 1.0 300.0 1.0 0.0 0.0
!      0.4 1.0 1.0 300.0 1.0 0.0 0.0
!      0.6 1.0 1.0 300.0 1.0 0.0 0.0
!      0.8 1.0 1.0 300.0 1.0 0.0 0.0
!      1.0 1.0 1.0 300.0 1.0 0.0 0.0
!     TABLE_DATA END
!     !- end of streams.in file
!     
! 5. Template for streams.in (modes=6,7,8,9)
!
!    !----------------------------------------------------------------------------!
!    !                                                                            !
!    !      Template for streams.in file (new format)                             !
!    !      (for modes 6, 7, 8, and 9)                                            !
!    !      HW, March, 2010                                                       !
!    !                                                                            !
!    !      Note: * any part beginning with "!" is a comment                      !
!    !            * empty lines are ignored                                       !
!    !            * keywords are not case-sensitive                               !
!    !                                                                            !
!    !----------------------------------------------------------------------------!
!    !                                                                   
!    ! specify MODECI: the mode of thermochemistry to be used in ISAT
!    !        = 1 or CONST_PROPERTY     - inert, constant-density flow
!    !        = 2 or MIXTURE_FRACTION   - mixture fraction formulation 
!    !        = 3 or PROGRESS_VARIABLE  - reaction progress variable formulation 
!    !        = 4                       - not supported yet
!    !        = 5                       - not supported yet
!    !        = 6 or DIRECT_INTEGRATION - CHEMKIN/direct integration
!    !        = 7 or ISAT_DI            - CHEMKIN/ISAT
!    !        = 8 or RCCE               - CHEMKIN/ISAT/RCCE
!    !        = 9 or ICE_PIC            - CHEMKIN/ISAT/ICE_PIC
!    !              
!    MODECI         ISAT_DI
!    !
!    !      Specify the thermodynamic state of each stream.
!    !      - Each stream begins with 'STREAM BEGIN' and ends with 'STREAM END'.
!    !      - Each stream should have a unique name, e.g., "FUEL", "PILOT", "COFLOW"
!    !      - The species can be specified in two units: [MOLE] and [MASS].
!    !      - [EQUIL] [EQUIL-H] [EQUIL-T] indicate to use the equilibrium state of the stream
!    !        [EQUIL] or [EQUIL-H]: constant enthalpy; [EQUIL-T]: constant temperature
!    !      - P  : pressure in atmosphere unit
!    !      - T  : temperature in unit [K] 
!    !
!    ! first stream: PILOT 
!    STREAM BEGIN PILOT [MOLE] [EQUIL-T] ! use the equilibrium state for the PILOT stream
!    	 P          0.993                ! pressure (atm)
!        T          1880                 ! temperature (K)
!        H2         0.0064               ! compositions in relative [MOLE] units
!        H          0.0025
!        O2         0.1699
!        OH         0.0165
!        H2O        0.5229
!        CO         0.0145
!        CO2        0.2495
!        N2         2.6209
!    STREAM END PILOT
!    ! second stream: COFLOW
!    STREAM BEGIN COFLOW [MASS]        ! use the stream as stated (not equilibrium)
!        P          0.993 
!        T          291
!    	 O2         21.                ! compositions in relative [MASS] units
!    	 N2         78.
!        H2O        1.00
!    STREAM END COFLOW
!    ! third stream: FUEL
!    STREAM BEGIN FUEL [MOLE]
!        P          0.993
!        T          294
!        O2         15.75
!        CH4        25.0
!        N2         59.25
!    STREAM END FUEL
!    !
!    !      Specify the represented species for RCCE (8) and ICE_PIC (9)
!    !      - Species names are enclosed between "REPRESENTED_SPECIES BEGIN" 
!    !        and "REPRESENTED_SPECIES END"
!    !      - Species names are separated by blanks or line breaks
!    !      - The input has no effect for MODECI<8 
!    !
!    REPRESENTED_SPECIES BEGIN
!    	CH4 O2 CO2 H2O 
!    	CO  H2 
!    	OH  H  N2
!    REPRESENTED_SPECIES END	
!    !- end of streams.in file
!
!END-OF-GUIDE
!
module streams_mod
    implicit none
!------------- private data begin ------------------------------------------------- 
    private
    integer,parameter        :: k_dp = kind(1.d0)   ! kind of double precision
    integer,parameter        :: MaxStrLen = 1000    ! max length of character string
    character(len=200)       :: CurrentLine=''      ! current line of streams file
    character(len=200)       :: ErrMsg=''           ! error message
    integer                  :: nStreamRead=-1      ! number of streams that have been read in
    character(len=16)        :: StreamNames         ! the name of streams
    allocatable              :: StreamNames(:)
    character(len=16)        :: Cnames              ! [6-9] names of full variables
    allocatable              :: Cnames(:)
    integer                  :: nC=-1               ! [6-9] number of full variables
    integer,parameter        :: imode = 1           ! modeci line 
    integer,parameter        :: istrm = 2           ! stream line
    integer,parameter        :: iresp = 3           ! represented_species line
    integer,parameter        :: idens = 4           ! density line
    integer,parameter        :: itabv = 5           ! table_variables line
    integer,parameter        :: itabd = 6           ! table_data line
    integer,parameter        :: itmsc = 7           ! timescale line  
    logical                  :: Initialized=.false. ! initialization status of this module
!------------- private data end --------------------------------------------------- 
    
!------------- public data begin --------------------------------------------------
    type,public :: StrmPublicType                   ! [modes]
        integer              :: iFormat=-1          ! [all] format of the streams file
        integer              :: ModeCI=-1           ! [all] ISAT-CK mode
        real(k_dp)           :: density=-1.d0       ! [1  ] density
        real(k_dp)           :: timescale=-1.d0     ! [3  ] timescale
        integer              :: nStream=-1          ! [2-9] number of streams
        integer              :: nC     =-1          ! [2-9] size of each Cvalue
        real(k_dp),pointer   :: Cvalue(:,:)=>NULL() ! [2-9] stream data (size=(nC,nStream))
        integer,pointer      :: kstrm(:)=>NULL()    ! [6-9] need to calculate equilibrium state of the stream if ==2 (const enth) or 3 (const temp)
        integer,pointer      :: kmole(:)=>NULL()    ! [2-9] units (=1:mole;=0:mass) of the returned stream compositions
        integer              :: nSp=-1              ! [2-9] number of species (represented species for mode=8,9)
        character*16,pointer :: Spnames(:)=>NULL()  ! [2-9] name of species (size=nSp) (represented species for mode=8,9)
        integer              :: nTable =-1          ! [2-3] size of table data (# rows)
        real(k_dp),pointer   :: Table(:,:)=>NULL()  ! [2-3] table data (size=(nTable,nSp+4))
    end type StrmPublicType
    type(StrmPublicType ),public  :: StrmDat
    
    public                   :: streams_init
    public                   :: streams_status
    public                   :: streams_free
    public                   :: streams_index
    public                   :: streams_name
    public                   :: streams_ErrMsg
    public                   :: streams_format
    public                   :: streams_modeci
    
! for testing use only ---------------
    public                   :: streams_nstream
    public                   :: streams_ntable
    public                   :: UpperCase    
    public                   :: CheckLine
    public                   :: StringRead
    public                   :: ModifyLine
    public                   :: LineType
!------------- public data end ----------------------------------------------------

contains
!------------------------------------------------------------------------------
subroutine streams_init(Unit,err,NAMES,NNAMES)             ! [public] initialize the streams module
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams file
    character*(*),intent(in) :: NAMES(:)           ! [optional] required for modes=2-9 only (names of full variables)
    integer,intent(in)       :: NNAMES             ! [optional] required for modes=2-9 only (number of full variables)
    optional                 :: NAMES,NNAMES
    logical,intent(out)      :: err                ! error
! - local variables
    logical                  :: L1,L2
    
! - check Unit    
    if(.not.CheckUnit(Unit)) then
        err = .true.
        return
    endif

! - checking    
    if(Initialized) then
        err = .true.
        CurrentLine = ''
        ErrMsg = 'streams module should be initialized once only (free it first) [streams_init]'
        return
    endif

! - input checking
    L1 = present(NAMES)
    L2 = present(NNAMES)
    if(L1.and.L2) then
        ! both present
        nC     = NNAMES
        allocate(Cnames(nC))
        Cnames = NAMES
    elseif(.not.(L1.or.L2)) then
        !  none present
        nC = -1
    else
        err = .true.
        CurrentLine = ''
        ErrMsg = 'Presenting only one of NAMES and NNAMES is not allowed for [streams_init]'
        return
    endif
    
! - format of streams file    
    StrmDat%iFormat = streams_format(Unit,err)
    if(err) return
    
! - read streams file
    if(StrmDat%iFormat==1) then
        ! old format
        !call streams_read_fmt1(Unit,err)
        err = .true.
        ErrMsg = 'reading old format has not been supported yet'
        return
    else
        ! new format
        call streams_read_fmt2(Unit,err)
    endif
    
    if(err) return
    
    Initialized = .true.    
    
end subroutine streams_init
!------------------------------------------------------------------------------
function   streams_status() result(init_status)            ! [public] return the initialization status of streams module
    implicit none
    logical                  :: init_status
    init_status = Initialized
end function streams_status
!------------------------------------------------------------------------------
subroutine streams_free()                                  ! [public] free the initialization
    implicit none
    StrmDat%iFormat=-1
    StrmDat%nStream=-1          
    StrmDat%nC     =-1          
    StrmDat%nTable =-1
    if(associated(StrmDat%Cvalue)) deallocate(StrmDat%Cvalue)
    if(associated(StrmDat%Table))  deallocate(StrmDat%Table)
    if(associated(StrmDat%kstrm))  deallocate(StrmDat%kstrm)
    if(associated(StrmDat%kmole))  deallocate(StrmDat%kmole)
    StrmDat%ModeCI=-1           
    StrmDat%density=-1.d0       
    StrmDat%timescale=-1.d0     
    StrmDat%nSp=-1              
    if(associated(StrmDat%Spnames)) deallocate(StrmDat%Spnames)
    Initialized = .false.
end subroutine streams_free
!------------------------------------------------------------------------------
function   streams_index(name,err) result(Iresult)         ! [public] return the stream index given stream name
    implicit none
    character*(*),intent(in) :: name
    logical,intent(out)      :: err
    integer                  :: Iresult
! - local variables

    Iresult = 0
    if(.not.initialized) then
        err = .true.
        ErrMsg = 'The streams are not initialized or the old format is used  [streams_index]'
        CurrentLine = ''
        return
    endif

    Iresult = StringMatch(name,StreamNames,nStreamRead)
    
end function streams_index
!------------------------------------------------------------------------------
function   streams_name(istream,err) result(Cresult)       ! [public] return the stream name given stream index
    implicit none
    integer,intent(in)       :: istream
    logical,intent(out)      :: err
    character*(16)           :: Cresult
! - local variables

    Cresult = ''
    if(.not.initialized) then
        err = .true.
        ErrMsg = 'The streams are not initialized or the old format is used [streams_name]'
        CurrentLine = ''
        return
    endif

    Cresult = StreamNames(istream)
    
end function streams_name
!------------------------------------------------------------------------------
subroutine streams_ErrMsg(line,error)                      ! [public] retrieve error message if err=.true.
    implicit none
    character*(*),intent(out):: line      ! line that contains the error
    character*(*),intent(out):: error     ! error message
    optional                 :: line,error
    
    if(present(line )) line = CurrentLine
    if(present(error)) error= ErrMsg
    
end subroutine streams_ErrMsg
!------------------------------------------------------------------------------
function   streams_format(Unit,err) result(Lformat)        ! [public] return the format of streams.in
! check the format of the streams.in file, and return 
!     1: for the old format 
!     2: for the new format (see streams.in_template for the new format)
! note: * this function does NOT check the correctness of the input in streams.in
!       * for the new format, more checking is done during streams_read
!
    implicit none
    integer,intent(in)       :: Unit     ! unit connecting to the streams.in file
    logical,intent(out)      :: err
    optional                 :: Unit
    integer                  :: Lformat  ! output the format
! - local variables
    character(len=500)       :: string
    logical                  :: eof,skip
    logical,save             :: isCalled=.false.
    integer,save             :: Lformat_saved
    
    if(isCalled.and. .not.present(Unit)) then
        Lformat = Lformat_saved
        err     = .false.
        return
    endif
    
    isCalled = .true.

    Lformat = 0
    
    if(.not.(present(Unit))) then
        ErrMsg   = 'Error: missing argument (Unit) when calling streams_format for the first time'
        isCalled = .false.
        err      = .true.
        return
    endif
    
    if(.not.CheckUnit(Unit)) then
        err = .true.
        return
    endif

    rewind(Unit)
    call GetOneLine(Unit,string,eof,err,skip)
    if(eof) then
        ErrMsg = 'empty streams.in file [streams_format]'
        err    = .true.
    endif
    if(err) return
    
    if(len_trim(string)==0) then
        ! The first line of streams.in is an empty line, so it can NOT be the old format
        Lformat = 2
        Lformat_saved = Lformat
        return
    else
        ! check the first non-blank character
        string = adjustl(string)
        if(ichar(string(1:1))>=ichar('0').and.ichar(string(1:1))<=ichar('9')) then
            ! The fisrt character is a number between 0-9, so it can NOT be the new format
            Lformat = 1
            Lformat_saved = Lformat
        else
            Lformat = 2
            Lformat_saved = Lformat
        endif
    endif
    
end function streams_format
!------------------------------------------------------------------------------
function   streams_modeci(Unit,err) result(modeci)         ! [public] read modeci from streams file
                        
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
    integer                  :: modeci
! - local variables
    integer                  :: i,I1
    character(len=MaxStrLen) :: line
    character(len=16)        :: string
    logical                  :: L1
    logical                  :: eof
    integer,parameter        :: nwk=10
    integer                  :: Iwk(nwk)
    real(k_dp)               :: Rwk(nwk)
    character(len=32)        :: Cwk(nwk)

    modeci = -1
    err    = .false.
    CurrentLine = ''
    
! - check unit
    if(.not.CheckUnit(Unit)) then
        err = .true.
        return
    endif
    
! - format of the streams file
    StrmDat%iFormat = streams_format(Unit,err)
    if(err) return
    
! - start to read
    rewind(Unit)
    if(StrmDat%iFormat==1) then
        ! old format
        read(Unit,*,err=100,end=101) modeci
        return
100     ErrMsg = 'error in reading modeci from the streams file [streams_modeci::1]'
        return
101     ErrMsg = 'can not read modeci from the streams file [streams_modeci::1]'
        return
    elseif(StrmDat%iFormat==2) then
        call GetInputLine(Unit,line,eof,err)
        do while(.not.(eof.or.err))
            i = LineType(line,Iwk,Rwk,Cwk,nwk,err)
            if(err) return
            if(i==imode) then
                modeci = Iwk(1)  
                exit
            endif
            call GetInputLine(Unit,line,eof,err)
        enddo
        if(err) return
        if(modeci==-1) then
            err = .true.
            ErrMsg = 'Can not read modeci from the streams file [streams_modeci::2]'
            return
        endif
    endif
        
end function streams_modeci
!------------------------------------------------------------------------------
!//////////////////////////////////////////////////////////////////////////
! The following are private subroutines and should be used internally only
!//////////////////////////////////////////////////////////////////////////
function   uppercase(str) result(str_uppercase)            ! [private] return uppercase of input string
    implicit none
    character*(*),intent(in) :: str
    character*(len(str))     :: str_uppercase
! - local variables
    integer                  :: i,m,ich,ia,iz,ishift
    
    ia = ichar('a')
    iz = ichar('z')
    ishift = ichar('A')-ia
    str_uppercase = str
    
    m = len(str)    
    do i=1,m
        ich = ichar(str(i:i))
        if(ich<=iz.and.ich>=ia) str_uppercase(i:i) = char(ich+ishift)
    enddo          
    
end function uppercase
!------------------------------------------------------------------------------
subroutine StringRead(Str,delimiter,StrArray,ns)           ! [private] read strings separated by delimiter
    implicit none
    character*(*),intent(in)    :: Str
    character*(*),intent(in)    :: delimiter
    character*(*),pointer       :: StrArray(:)
    integer,intent(out)         :: ns
! - local variables
    character(len=len(Str))     :: tline
    integer                     :: nvalue,nvalue_max
    character(len=len(StrArray)),pointer:: sarray(:),sarray_bak(:)
    integer                     :: len_a,len_d,i

    ! deallocate StrArray
    if(associated(StrArray)) deallocate(StrArray)
    
    nvalue_max = 100
    allocate(sarray(nvalue_max))
    
    ! input checking
    if(len_trim(Str)==0) then
        nvalue = 0
        return
    endif
    if(len(delimiter)==0) then
        sarray(1) = trim(adjustl(Str))
        nvalue    = 1
    else        
        tline = ''
        tline = adjustl(Str)
        len_a = len(tline)
        len_d = len(delimiter)
        nvalue= 0
        
        i     = index(tline,delimiter)
        do while(i>0.and.len_trim(tline)>0) 
            nvalue = nvalue+1
            if(nvalue>nvalue_max) then
                ! extend the space for sarray
                if(associated(sarray_bak)) deallocate(sarray_bak)
                allocate(sarray_bak(nvalue_max*2))
                sarray_bak = ''
                sarray_bak(1:nvalue_max) = sarray
                deallocate(sarray)
                sarray => sarray_bak
                nvalue_max = nvalue_max*2
                sarray_bak => NULL()
            endif
            sarray(nvalue) = tline(1:i-1)
            tline(1:i+len_d-1) = ''
            tline = adjustl(tline)
            i     = index(trim(tline),delimiter)
        enddo
        if(len_trim(tline)/=0) then
            nvalue = nvalue+1
            sarray(nvalue) = trim(adjustl(tline))
        endif
    endif
    
    ns = nvalue
    allocate(StrArray(ns))
    StrArray(1:ns) = sarray(1:ns)
    
    deallocate(sarray)        

end subroutine StringRead
!------------------------------------------------------------------------------
function   StringMatch(str,StrArray,ns) result(imatch)     ! [private] find a match of str in StrArray (not case-sensitive)
    implicit none
    character*(*),intent(in)    :: Str
    character*(*),intent(in)    :: StrArray(ns)
    integer,intent(in)          :: ns
    integer                     :: imatch
! - local variables
    integer                     :: i

    imatch = 0
    do i=1,ns
        if(uppercase(Str)==uppercase(StrArray(i))) then
            imatch = i
            exit
        endif
    enddo

end function StringMatch
!------------------------------------------------------------------------------
function   CheckUnit(Unit) result(Lresult)                 ! [private] check Unit connecting to a valid file or not
    implicit none
    integer,intent(in)       :: Unit
    logical                  :: Lresult

    Inquire(UNIT=Unit,OPENED=Lresult)
    if(Lresult) ErrMsg = 'Unit is NOT connecting to a streams.in file [CheckUnit]'
    
end function CheckUnit
!------------------------------------------------------------------------------
function   CheckLine(line,cmark) result(Iresult)           ! [private] check for comment line or empty line
    implicit none
    character*(*),intent(in) :: line     ! input line
    character,intent(in)     :: cmark    ! character indicating comment line
    integer                  :: Iresult  ! 1: comment line; 2: empty line; 0: none of the above
! - local variables
    character(len=len_trim(line)):: line2
    
    Iresult = 0
    
    if(len_trim(line)<=0) then
        Iresult = 2
    else
        line2 = adjustl(line)
        if(line2(1:1)==cmark) Iresult = 1
    endif
    
end function CheckLine
!------------------------------------------------------------------------------
subroutine ModifyLine(line,line_new,icase)                 ! [private] modify the line (replace Tab and remove comment)
    implicit none
    character*(*),intent(inout) :: line
    character*(*),intent(inout) :: line_new
    integer,intent(in)       :: icase  ! 1: replace Tab with blank; 
                                       ! 2: remove comment line starting with '!';
                                       ! 0: all of the above
! - local variables
    integer                  :: i,m                                       
    character                :: CTab
    
    line_new = line
    
    CTab = char(9)
    if(icase==0.or.icase==1) then
        m = len_trim(line_new)
        do i=1,m
            if(line_new(i:i)==CTab) line_new(i:i) = ' '
        enddo
    endif
    
    if(icase==0.or.icase==2) then
        i = index(line_new,'!')
        if(i>0) line_new = line_new(1:i-1)        
    endif
        
end subroutine ModifyLine
!------------------------------------------------------------------------------
subroutine GetOneLine(Unit,line,eof,err,skip)              ! [private] get one line
    implicit none
    integer,intent(in)       :: Unit
    character*(*),intent(out):: line
    logical,intent(out)      :: eof,err,skip
    
    eof = .false.
    err = .false.
    skip= .true.
    read(Unit,'(a)',ERR=100,END=200) CurrentLine
    call ModifyLine(CurrentLine,line,0)
    if(CheckLine(line,'!')==0) skip=.false.
    
    return
    
100 err = .true.
    return
200 eof = .true.
    return
    
end subroutine GetOneline
!------------------------------------------------------------------------------
subroutine GetInputLine(Unit,line,eof,err)                 ! [private] get next input line (skip comments and empty lines)
    implicit none
    integer,intent(in)       :: Unit
    character*(*),intent(out):: line
    logical,intent(out)      :: eof,err
! - local variables    
    logical                  :: skip
    
    skip = .true.
    do while(skip)
        call GetOneLine(Unit,line,eof,err,skip)
        if(eof.or.err) return
    enddo
    
end subroutine GetInputLine
!------------------------------------------------------------------------------
function   LineType(line,Iwk,Rwk,Cwk,nwk,err) result(Itype)! [private] return the type of the line
    implicit none
    character*(*),intent(in) :: line  
    integer,intent(in)       :: nwk        ! size of workspace
    integer,intent(inout)    :: Iwk(nwk)   ! integer workspace
    real(k_dp),intent(inout) :: Rwk(nwk)   ! real workspace
    character*(*),intent(inout):: Cwk(nwk) ! character workspace
    logical,intent(out)      :: err
    integer                  :: Itype      
!IType= 1: 'MODECI ....'
!       2: 'STREAM BEGIN ....'
!      -2: 'STREAM END ...'                                     
!       3: 'REPRESENTED_SPECIES BEGIN ...'                                     
!      -3: 'REPRESENTED_SPECIES END ...'                                     
!       4: 'DENSITY 1.0'                                     
!       5: 'TABLE_VARIABLES BEGIN ...'                                     
!      -5: 'TABLE_VARIABLES END'                                     
!       6: 'TABLE_DATA BEGIN'                                     
!      -6: 'TABLE_DATA END'                                     
!       7: 'TIMESCALE 0.3'       
!       0: others                              
! - local variables
    integer                  :: ns,i,j
    character(len=32),pointer:: StrArray(:)    
    integer,parameter        :: nkeys = 7
    character(len=32)        :: keys(nkeys)
    integer,parameter        :: nmode = 9              ! number of supported modeci
    character(len=32)        :: cmode(nmode)           ! names of modeci
    data keys /'MODECI','STREAM','REPRESENTED_SPECIES','DENSITY','TABLE_VARIABLES',&
               'TABLE_DATA','TIMESCALE'/
    data cmode /'CONST_PROPERTY','MIXTURE_FRACTION','PROGRESS_VARIABLE' , &
                'NULL:MODECI=4' ,'NULL:MODECI=5'   ,'DIRECT_INTEGRATION', &
                'ISAT_DI'       ,'RCCE'            ,'ICE_PIC'/
    
    Itype   = 0
    err     = .false.
    
    ! separate keywords
    call StringRead(line,' ',StrArray,ns)
    
    if(ns<2) return   ! not the line we want to find
    
    ! match the keywords
    Itype = StringMatch(StrArray(1),keys,nkeys)
    
    ! error checking
    if(Itype<1) then
        Itype = 0
        return
    endif
    
    ! check workspace (require nwk>=8 for safety)
    if(nwk<8) then
        err         = .true.
        CurrentLine = ''
        ErrMsg      = 'workspace must be at least 8 [LineType]'
        return
    endif
    Iwk = -10000000
    Rwk = -1.D-20
    Cwk = ''
    
    ! treat different types of lines separately
    select case (Itype)
        !_______________________________________________________________________________    
        !//////////////// MODECI line: modeci=Iwk(1)
        case (imode)    
            if(len_trim(StrArray(2))>0) then
                if(ichar(StrArray(2)(1:1))<ichar('0').or.ichar(StrArray(2)(1:1))>ichar('9')) then
                    ! mode name is detected
                    Iwk(1) = StringMatch(Uppercase(StrArray(2)),cmode,nmode)
                    if(Iwk(1)<1) then
                        err = .true.
                        ErrMsg = 'do not recognize the MODECI name [LineType]'
                    endif
                else
                    ! mode number is detected
                    read(StrArray(2),*,ERR=102) Iwk(1)
                endif
            else
                err    = .true.
                ErrMsg = 'No MODECI value is specified on the MODECI line [LineType]'
            endif
            return
        !_______________________________________________________________________________    
        !//////////////// STREAM BEGIN/END line: Cwk(1)=Stream_Name, Data Unit=Iwk(1), Equil=Iwk(2)
        case (istrm)   
            if(ns<3) then
                err         = .true.
                CurrentLine = trim(line)
                ErrMsg      = 'missing stream name in the STREAM BEGIN/END line [LineType]'
                return
            endif
            Cwk(1) = StrArray(3)
            if(uppercase(StrArray(2))=='BEGIN') then
                ! find optional input '[MOLE]','[MASS]','[EQUIL]','[EQUIL-H]','[EQUIL-T]'
                Iwk(1) = 1     ! mole unit in default
                Iwk(2) = 1     ! not equilibrium in default
                do j=4,ns
                    if(uppercase(StrArray(j))=='[MOLE]') then
                        Iwk(1) = 1
                    elseif(uppercase(StrArray(j))=='[MASS]') then
                        Iwk(1) = 0
                    elseif(uppercase(StrArray(j))=='[EQUIL]'.or.uppercase(StrArray(j))=='[EQUIL-H]') then
                        Iwk(2) = 2
                    elseif(uppercase(StrArray(j))=='[EQUIL-T]') then
                        Iwk(2) = 3
                    else
                        CurrentLine = trim(line)
                        ErrMsg      = 'Do not recognize keyword '//trim(StrArray(j))//' on the line [LineType]'
                        err         = .true.
                        return
                    endif
                enddo
            elseif(uppercase(StrArray(2))=='END') then
                Itype = -Itype
            else
                err    = .true.
                ErrMsg = 'Do not recognize keyword '//trim(StrArray(2))//' on the line [LineType]'
                CurrentLine = line
            endif
            return
        !_______________________________________________________________________________    
        !//////////////// REPRESENTED_SPECIES/TABLE_VARIABLES BEGIN/END line
        case (iresp, itabv)   
            if(uppercase(StrArray(2))=='BEGIN') then
                ! do nothing
            elseif(uppercase(StrArray(2))=='END') then
                Itype = -Itype
            else
                err    = .true.
                ErrMsg = 'Do not recognize keyword '//trim(StrArray(2))//' on the line [LineType]'
                CurrentLine = line
            endif         
            return   
        !_______________________________________________________________________________    
        !//////////////// DENSITY/TIMESCALE line: DENSITY/TIMESCALE=Rwk(1)
        case (idens,itmsc)    
            if(len_trim(StrArray(2))>0) then
                read(StrArray(2),*,ERR=103) Rwk(1)
            else
                err    = .true.
                ErrMsg = 'No DENSITY/TIMESCALE value is specified on the line [LineType]'
            endif
            return
        !_______________________________________________________________________________    
        !//////////////// TABLE_DATA BEGIN/END line: Data Unit=Iwk(1)
        case (itabd)   
            if(uppercase(StrArray(2))=='BEGIN') then
                ! find optional input '[MOLE]','[MASS]'
                Iwk(1) = 1     ! mole unit in default
                do j=3,ns
                    if(uppercase(StrArray(j))=='[MOLE]') then
                        Iwk(1) = 1
                    elseif(uppercase(StrArray(j))=='[MASS]') then
                        Iwk(1) = 0
                    else
                        CurrentLine = trim(line)
                        ErrMsg      = 'Do not recognize keyword '//trim(StrArray(j))//' on the line [LineType]'
                        err         = .true.
                        return
                    endif
                enddo
            elseif(uppercase(StrArray(2))=='END') then
                Itype = -Itype
            else
                err    = .true.
                ErrMsg = 'Do not recognize keyword '//trim(StrArray(2))//' on the line [LineType]'
                CurrentLine = line
            endif
            return
        !_______________________________________________________________________________    
        !//////////////// error
        case default
            CurrentLine = ''
            ErrMsg      = 'unexpected value of Itype [LineType]'
            err         = .true.
            return
    end select
    return

102 err = .true.
    ErrMsg = 'Error in reading the MODECI line [LineType]'
103 err = .true.
    ErrMsg = 'Error in reading the DENSITY/TIMESCALE line [LineType]'    
    
end function LineType
!------------------------------------------------------------------------------
function   LineTypeQ(line,err) result(Itype)               ! [private] simplified version of LineType
    implicit none
    character*(*),intent(in) :: line  
    logical,intent(out)      :: err
    integer                  :: Itype      
! - local variables
    integer,parameter        :: nwk = 10
    integer                  :: Iwk(nwk)   
    real(k_dp)               :: Rwk(nwk)   
    character(len=32)        :: Cwk(nwk) 

    Itype = LineType(line,Iwk,Rwk,Cwk,nwk,err)
    
end function LineTypeQ
!------------------------------------------------------------------------------
function   streams_nstream(Unit,err) result(nstr)          ! [private] return number of streams
    implicit none
    integer,intent(in)       :: Unit     ! unit connecting to the streams.in file
    integer                  :: nstr     ! number of streams
    logical,intent(out)      :: err      ! .true. if an error occurs
! - local variables
    character(len=1000)      :: line
    character(len=16)        :: string
    logical                  :: eof,L1
    integer                  :: I1

! - check file unit
    if(.not.CheckUnit(Unit)) then
        err = .true.
        return
    endif

! - find number of streams ('STREAM BEGIN ...' for modes=2-3, 6-9)
    rewind(Unit)
    call GetInputLine(Unit,line,eof,err)
    nstr = 0
    do while(.not.(eof.or.err))
        if(LineTypeQ(line,err)==istrm) nstr = nstr+1
        if(err) exit
        call GetInputLine(Unit,line,eof,err)
    enddo
    
    if(err) then
        nstr = -1
        return
    endif
    
end function streams_nstream
!------------------------------------------------------------------------------
function   streams_ntable(Unit,err) result(nstr)           ! [private] return table size (number of rows)
    implicit none
    integer,intent(in)       :: Unit     ! unit connecting to the streams.in file
    integer                  :: nstr     ! number of streams
    logical,intent(out)      :: err      ! .true. if an error occurs
! - local variables
    character(len=1000)      :: line
    character(len=16)        :: string
    logical                  :: eof,L1
    integer                  :: I1

! - check file unit
    if(.not.CheckUnit(Unit)) then
        err = .true.
        return
    endif

! - test the 'TABLE_DATA BEGIN'
    rewind(Unit)
    call GetInputLine(Unit,line,eof,err)
    nstr = 0
    do while(.not.(eof.or.err))
        if(LineTypeQ(line,err)==itabd) exit
        if(err) exit
        call GetInputLine(Unit,line,eof,err)
    enddo
    if(err) then
        nstr = -1
        return
    endif

! - read until 'TABLE_DATA END'
!   note: no non-data lines are allowd in the TABLE_DATA section
    call GetInputLine(Unit,line,eof,err)
    do while(.not.(eof.or.err))
        if(LineTypeQ(line,err)==-itabd) exit
        if(err) exit
        nstr = nstr + 1
        call GetInputLine(Unit,line,eof,err)
    enddo
    
    if(err) then
        nstr = -1
        return
    endif
    
end function streams_ntable
!------------------------------------------------------------------------------
subroutine ReadStream(Unit,nC,Cnames,Cvalue, &             ! [private] read current stream
                      CurrentStrmName,err)
    implicit none
    integer,intent(in)       :: Unit
    integer,intent(in)       :: nC                 ! number of compositions
    character*(*),intent(in) :: Cnames(nC)         ! names of compositions
    real(8),intent(inout)    :: Cvalue(nC)         ! composition in relative mass or mole units
                                                   ! for compositions not specified in streams.in, their values are unchanged
    character*(*),intent(in) :: CurrentStrmName    ! current stream name
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
! - local variables
    integer                  :: i,nread
    character(len=MaxStrLen) :: line
    character(len=16)        :: string
    logical                  :: eof    
    integer                  :: ns
    character(len=16),pointer:: StrArray(:)    
    integer,parameter        :: nkeys=2            ! number of required key input (P and T for modes 6-9)
                                                   ! not count the species and mixture fraction/progress variable
    character(len=16)        :: keys(nkeys)
    logical                  :: FindKeys(nkeys)
    data keys /'P','T'/

    FindKeys    = .false.
    nread       = 0

    call GetInputLine(Unit,line,eof,err)
    do while(.not.(eof.or.err))
        if(LineTypeQ(line,err)==-Istrm) then
            if(err) return
            i = index(line,CurrentStrmName)
            if(i<=0) then
                err = .true.
                CurrentLine = line
                ErrMsg      = 'Stream name does not match '//CurrentStrmName
                return
            endif 
            exit
        else
            if(err) exit
            ! read key and value
            call StringRead(line,' ',StrArray,ns)
            if(ns/=2) then
                err    = .true.
                ErrMsg = 'Do not recognize the line (KEY  VALUE) [ReadStream]'
                exit
            endif
            i = StringMatch(StrArray(1),Cnames,nC)
            if(i<1) then
                err    = .true.
                ErrMsg = 'Do not recognize the key '//trim(StrArray(1))//' [ReadStream]'
                exit
            endif
            if(len_trim(StrArray(2))>0) then
                read(StrArray(2),*,ERR=100) Cvalue(i)
                nread = nread+1
            else
                err = .true.
                ErrMsg = 'Invalid input line in streams.in [ReadStream]'
                exit
            endif
            ! required keys
            i = StringMatch(StrArray(1),keys,nkeys)
            if(i>0) FindKeys(i)=.true.
        endif
        call GetInputLine(Unit,line,eof,err)
    enddo
    
    if(err) return
    
    ! check required keys
    ! for modes = 2,3, P and T are the required keys
    if(StrmDat%ModeCI==2.or.StrmDat%ModeCI==3) then
        ! set FindKeys=.true. without checking for modes=2,3
        i = StringMatch('T',keys,nkeys) 
        if(i>0) FindKeys(i)=.true.
        nread = nread + 1   ! T
        i = StringMatch('P',keys,nkeys) 
        if(i>0) FindKeys(i)=.true.
        nread = nread + 1   ! P
    endif
    ! for modes = 6-9, P, T is the required key
    do i=1,nkeys
        if(.not.findkeys(i)) then
            err = .true.
            ErrMsg = 'The key '//trim(keys(i))//' is missing in stream '//trim(CurrentStrmName)//' [ReadStream]'
            return
        endif
    enddo
    
    if(nread==nkeys) then
        err = .true.
        ErrMsg = 'No species or "f" exist in input stream '//trim(CurrentStrmName)//' [ReadStream]'
        return
    endif
    
    return
    
100 ErrMsg = 'Error in reading the value for the key '//trim(StrArray(1))//' [ReadStream]'
    err = .true.
                
end subroutine ReadStream 
!------------------------------------------------------------------------------
subroutine ReadSpeciesNames(Unit,nSp,Spnames,err)          ! [private] read species names
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    integer,intent(out)      :: nSp                ! number of species
    character*(*),intent(out):: Spnames            ! name of species (allocation is done inside)
    pointer                  :: Spnames(:)         ! size=nSp when return
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
! - local variables
    integer,parameter        :: MaxL = MaxStrLen*10
    integer                  :: i
    character(len=MaxL)      :: BigLine
    character(len=MaxStrLen) :: line
    character(len=16)        :: string
    logical                  :: eof    

    BigLine = ''

    call GetInputLine(Unit,line,eof,err)
    do while(.not.(eof.or.err))
        if(LineTypeQ(line,err)==-Iresp.or.LineTypeQ(line,err)==-Itabv) then
            exit
        else
            if(err) exit
            if(len_trim(BigLine)+len_trim(line)+1>MaxL) then
                err = .true.
                ErrMsg = 'Too many species to read (increase MaxL) [ReadSpeciesNames]'
                exit
            endif
            BigLine = trim(BigLine)//' '//trim(line)
        endif
        call GetInputLine(Unit,line,eof,err)
    enddo
    
    if(err) return

    call StringRead(BigLine,' ',Spnames,nSp)

end subroutine ReadSpeciesNames 
!------------------------------------------------------------------------------
subroutine ReadTableData(Unit,value,n1,n2,err)             ! [private] read 2D table data
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    integer,intent(in)       :: n1,n2
    real(k_dp),intent(out)   :: value(n1,n2)
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
! - local variables
    integer                  :: i

    err = .false.
    
    do i=1,n2
        read(Unit,*,err=104,end=105) value(:,i)
    enddo
    return
    
104 err = .true.
    CurrentLine = ''
    ErrMsg = 'error in reading table data'    
    return
105 err = .true.
    CurrentLine = ''
    ErrMsg = 'not enough table data to read'    
    return

end subroutine ReadTableData
!------------------------------------------------------------------------------
subroutine streams_read_fmt2(Unit,err)                     ! [private] read new format of streams file
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    

    err = .false.

! - double check the format of the streams file
    if(StrmDat%iFormat<0) then
        StrmDat%iFormat = streams_format(Unit,err)
        if(err) return
    endif
    if(StrmDat%iFormat/=2) then
        err = .true.
        CurrentLine = ''
        ErrMsg = 'streams_read_fmt2 is not for old streams format'
        return
    endif

! - find modeci
    if(StrmDat%ModeCI<0) then
        StrmDat%ModeCI = streams_modeci(Unit,err)
        if(err) return
    endif
    
! - find number of streams
    StrmDat%nStream = streams_nstream(Unit,err)
    if(err) return
    if(StrmDat%nStream>0) then
        if(associated(StrmDat%kstrm)) deallocate(StrmDat%kstrm)
        if(associated(StrmDat%kmole)) deallocate(StrmDat%kmole)
        allocate(StrmDat%kstrm(StrmDat%nStream))
        allocate(StrmDat%kmole(StrmDat%nStream))
        StrmDat%kstrm = -1
        StrmDat%kmole = -1
        if(nC>0) then
            StrmDat%nC = nC
            if(associated(StrmDat%Cvalue)) deallocate(StrmDat%Cvalue)
            allocate(StrmDat%Cvalue(StrmDat%nC,StrmDat%nstream))
            StrmDat%Cvalue = 0.d0
        endif
    endif
    
! - find table size
    StrmDat%nTable = streams_ntable(Unit,err)    
    
! - treat different modes separately
    select case(StrmDat%ModeCI)
        case (1)
            call streams_read_fmt2_1(Unit,err)
            if(err) return
        case (2,3)    
            call streams_read_fmt2_23(Unit,err)
            if(err) return
        case (6,7,8,9)
            call streams_read_fmt2_6789(Unit,err)
            if(err) return
        case default
            err = .true.
            CurrentLine = ''
            ErrMsg = 'ModeCI is not supported'
            return
    end select
 
end subroutine streams_read_fmt2
!------------------------------------------------------------------------------
subroutine streams_read_fmt2_1(Unit,err)                   ! [private] read new format of streams file for mode=1
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
! - local variables
    integer                  :: i,I1
    character(len=MaxStrLen) :: line
    character(len=16)        :: string
    logical                  :: L1
    logical                  :: eof
    integer,parameter        :: nwk=10
    integer                  :: Iwk(nwk)
    real(k_dp)               :: Rwk(nwk)
    character(len=32)        :: Cwk(nwk)

    err = .false.
    
    rewind(Unit)
    call GetInputLine(Unit,line,eof,err)
    do while(.not.(eof.or.err))
        if(LineType(line,Iwk,Rwk,Cwk,nwk,err)==imode) then
            ! read MODECI
            if(err) return
            if(Iwk(1)/=StrmDat%MODECI) then
                CurrentLine = line
                ErrMsg      = 'inconsistent MODECI [streams_read_fmt2_1]'
                err         = .true.
                return
            endif
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==idens) then
            ! read density
            if(err) return
            StrmDat%density = Rwk(1)
        else
            err    = .true.
            ErrMsg = 'Do not recognize the input line in streams.in [streams_read_fmt2_1]'
            return
        endif
        call GetInputLine(Unit,line,eof,err)
    enddo

    if(err) return

    if(StrmDat%density==-1.d0) then
        err = .true.
        CurrentLine = ''
        ErrMsg = 'invalid or no input for density [streams_read_fmt2_1]'
        return
    endif

end subroutine streams_read_fmt2_1
!------------------------------------------------------------------------------
subroutine streams_read_fmt2_23(Unit,err)                  ! [private] read new format of streams file for mode=2,3
                        
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
! - local variables
    integer                  :: i,I1
    character(len=MaxStrLen) :: line
    character(len=16)        :: string
    logical                  :: L1
    logical                  :: eof
    integer,parameter        :: nwk=10
    integer                  :: Iwk(nwk)
    real(k_dp)               :: Rwk(nwk)
    character(len=32)        :: Cwk(nwk)

    err = .false.

    nStreamRead = 0
    if(allocated(StreamNames)) deallocate(StreamNames)
    if(StrmDat%nstream<1) then
        CurrentLine = ''
        ErrMsg      = 'Stream data is required for modes 2 or 3 [streams_read_fmt2_23]'
        err         = .true.
        return
    endif
    allocate(StreamNames(StrmDat%nstream))
    StreamNames = ''

    if(StrmDat%nC<1) then
        CurrentLine = ''
        ErrMsg      = 'number of variables is not correct [streams_read_fmt2_23]'
        err         = .true.
        return
    endif
    
    rewind(Unit)
    call GetInputLine(Unit,line,eof,err)
    do while(.not.(eof.or.err))
        if(LineType(line,Iwk,Rwk,Cwk,nwk,err)==imode) then
            ! read MODECI
            if(err) return
            if(Iwk(1)/=StrmDat%MODECI) then
                CurrentLine = line
                ErrMsg      = 'inconsistent MODECI [streams_read_fmt2_1]'
                err         = .true.
                return
            endif
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==itmsc) then
            ! read timescale (only for mode 3)
            if(err) return
            StrmDat%timescale = Rwk(1)
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==itabv) then
            ! read table variables (species names)
            if(err) return
            call ReadSpeciesNames(Unit,StrmDat%nSp,StrmDat%Spnames,err)
            if(err) return
            if(associated(StrmDat%Table)) deallocate(StrmDat%Table)
            allocate(StrmDat%Table(StrmDat%nSp+4,StrmDat%nTable))
            StrmDat%Table = 0.d0
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==itabd) then
            ! read table data
            if(err) return
            StrmDat%kmole = Iwk(1)
            call ReadTableData(Unit,StrmDat%Table,StrmDat%nSp+4,StrmDat%nTable,err)
            if(err) return
            call GetInputLine(Unit,line,eof,err)   ! get "TABLE_DATA END" line
            if(err) return
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==istrm) then
            ! streams
            nStreamRead = nStreamRead + 1
            StrmDat%Kmole(nStreamRead)  = Iwk(1)
            StrmDat%kstrm(nStreamRead)  = Iwk(2)
            StreamNames(nStreamRead) = trim(Cwk(1))
            call ReadStream(Unit,nC,Cnames,StrmDat%Cvalue(:,nStreamRead),StreamNames(nStreamRead),err)
            if(err) exit
        else
            err    = .true.
            ErrMsg = 'Do not recognize the input line in streams.in [streams_read_fmt2_23]'
            return
        endif
        call GetInputLine(Unit,line,eof,err)
    enddo

    if(err) return

    if(StrmDat%nC<=0.or.StrmDat%nStream<=0) then
        err = .true.
        CurrentLine = ''
        ErrMsg = 'Missing Table_Variables or Table_Data [streams_read_fmt2_23]'
        return
    endif
   
end subroutine streams_read_fmt2_23
!------------------------------------------------------------------------------
subroutine streams_read_fmt2_6789(Unit,err)                ! [private] read new format of streams file for mode=6,7,8,9
                        
    implicit none
    integer,intent(in)       :: Unit               ! unit connecting the streams.in file
    logical,intent(out)      :: err                ! an error occured during the call if .true. is returned    
! - local variables
    integer                  :: i,I1
    character(len=MaxStrLen) :: line
    character(len=16)        :: string
    logical                  :: L1
    logical                  :: eof
    integer,parameter        :: nwk=10
    integer                  :: Iwk(nwk)
    real(k_dp)               :: Rwk(nwk)
    character(len=32)        :: Cwk(nwk)

    err = .false.
    
    StrmDat%nSp = 0
    nStreamRead = 0
    if(allocated(StreamNames)) deallocate(StreamNames)
    if(StrmDat%nstream<1) then
        CurrentLine = ''
        ErrMsg      = 'Stream data is required for mode>6 [streams_read_fmt2_6789]'
        err         = .true.
        return
    endif
    allocate(StreamNames(StrmDat%nstream))
    StreamNames = ''

    if(StrmDat%nC<1) then
        CurrentLine = ''
        ErrMsg      = 'number of variables is not correct [streams_read_fmt2_6789]'
        err         = .true.
        return
    endif
    
    rewind(Unit)
    call GetInputLine(Unit,line,eof,err)
    do while(.not.(eof.or.err))
        if(LineType(line,Iwk,Rwk,Cwk,nwk,err)==imode) then
            ! read MODECI
            if(err) return
            if(Iwk(1)/=StrmDat%MODECI) then
                CurrentLine = line
                ErrMsg      = 'inconsistent MODECI [streams_read_fmt2_6789]'
                err         = .true.
                return
            endif
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==istrm) then
            ! streams
            nStreamRead = nStreamRead + 1
            StrmDat%Kmole(nStreamRead)  = Iwk(1)
            StrmDat%kstrm(nStreamRead)  = Iwk(2)
            StreamNames(nStreamRead) = trim(Cwk(1))
            call ReadStream(Unit,nC,Cnames,StrmDat%Cvalue(:,nStreamRead),StreamNames(nStreamRead),err)
            if(err) exit
        elseif(LineType(line,Iwk,Rwk,Cwk,nwk,err)==iresp) then
            call ReadSpeciesNames(Unit,StrmDat%nSp,StrmDat%Spnames,err)
            if(err) exit
            if(StrmDat%nSp>0) then
                ! check all represented species in Cnames
                do i=1,StrmDat%nSp
                    I1 = StringMatch(StrmDat%Spnames(i),Cnames,nC)
                    if(I1<1) then
                        err = .true.
                        ErrMsg = 'Represented species '//trim(StrmDat%Spnames(i))//' does not exist [streams_read_fmt2_6789]'
                        exit
                    else
                        ! use those in Cnames for RSnames to match the case 
                        StrmDat%Spnames(i) = Cnames(I1)
                    endif
                enddo
                if(err) exit
            endif
        else
            err    = .true.
            ErrMsg = 'Do not recognize the input line in streams.in [streams_read_fmt2_6789]'
            return
        endif
        call GetInputLine(Unit,line,eof,err)
    enddo

    if(err) return

    if(nStreamRead==0) then
        err = .true.
        ErrMsg = 'Can not find input stream compositions in streams.in [streams_read_fmt2_6789]'
        return
    endif
   
end subroutine streams_read_fmt2_6789
!------------------------------------------------------------------------------
end module streams_mod