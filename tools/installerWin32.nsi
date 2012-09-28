;---------------------------------------------------------------------------------------------
; installeurSofa.nsi
;
; Fait par Micha�l ADAM, le 12/08/08
;
; Ce script permet de configurer l'installation des fichiers et les options qui vont avec.
; Pour modifier ce fichier, installez NSIS et reportez vous � la doc et aux exemples.
;
;---------------------------------------------------------------------------------------------


;version de Sofa  =========> a mettre a jour!!!
!define VERSION "_1.0.3"


; Nom de l'installeur
Name "Sofa"

; Fichier d'installation executable
OutFile "Sofa${VERSION}-Setup.exe"

; r�pertoire d'installation par d�faut
InstallDir $PROGRAMFILES\Sofa${VERSION}
FileErrorText "An error occured during the copy of the files. $\n\
You may not have the permission to write on this location.$\n\
Contact your administrator for more details."

; cl� de registre pour se souvenir o� le logisiel a �t� install� (pour la d�sinstallation)
InstallDirRegKey HKLM "Software\Sofa" "Install_Dir"

Var modeler_installed
Var sofaCUDA_installed
Var doc_installed
Var examples_installed


;--------------------------------
; listes des pages utilis�es

Page license

Page components

Page directory

Page instfiles

UninstPage uninstConfirm
UninstPage instfiles


;--------------------------------
;----------- sp�cifie le fichier de licence � afficher
LicenseData ../LICENCE.txt
  
  InstType "Typical"
  InstType "Full"
  InstType "Base"  

;--------------------------------
;----------- Tous les trucs � installer

;---------
SectionGroup "Application"

  ;--------- executable et dlls
  Section "runSofa"
  ;----- runSofa est pr�sent dans les trois options d'installation (typical, full, base)
  SectionIn RO 1 2 3
  
    CreateDirectory $INSTDIR\bin
    SetOutPath $INSTDIR\bin
    File  "..\bin\runSofa.exe"
    File "..\bin\*.dll"

	;--- Example of the Liver
    CreateDirectory $INSTDIR\examples\Demos
    CreateDirectory $INSTDIR\examples\Simulation
    CreateDirectory $INSTDIR\share\BehaviorModels
    CreateDirectory $INSTDIR\share\config
    
    SetOutPath $INSTDIR\Share\config
    File /r "..\share\config\*" 
    CreateDirectory $INSTDIR\share\materials
    CreateDirectory $INSTDIR\share\mesh
    CreateDirectory $INSTDIR\share\screenshots
    
    CreateDirectory $INSTDIR\share\textures    
    SetOutPath $INSTDIR\Share\textures    
    File "..\share\textures\media*" 
    File "..\share\textures\SOFA_logo.bmp"
    
    SetOutPath $INSTDIR\share\mesh
    File "..\share\mesh\liver.msh" 
    File "..\share\mesh\liver.mtl" 
    File "..\share\mesh\liver-smooth.obj" 
    File "..\share\mesh\liver.sph" 
    
    
    SetOutPath $INSTDIR\examples\Demos
    File "..\examples\Demos\liver.scn"

    ;---------- �criture des infos dans les cl�s de registre
    ;--- dabord on v�rifie les droit de l'utilisateur : si admin->install pour tout le monde si simple user->install que pour lui
    Pop $0
    UserInfo::GetAccountType
    Pop $1
    UserInfo::GetOriginalAccountType
    Pop $2
    StrCmp $1 "Admin" 0 +3
      SetShellVarContext all
      Goto done
    StrCmp $1 "Power" 0 +3
      SetShellVarContext all
      Goto done
    StrCmp $1 "User" 0 +3
      SetShellVarContext current
      Goto done
    StrCmp $1 "Guest" 0 +3
      SetShellVarContext current
      Goto done
    MessageBox MB_OK "Unknown error"
    Goto done

    done:

    WriteRegStr HKLM SOFTWARE\Sofa "Install_Dir" "$INSTDIR"
  
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Sofa${VERSION}" "DisplayName" "Sofa"
    WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Sofa${VERSION}" "UninstallString" '"$INSTDIR\uninstall.exe"'
    WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Sofa${VERSION}" "NoModify" 1
    WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Sofa${VERSION}" "NoRepair" 1
    WriteUninstaller "uninstall.exe"

    ;--- essai pour l'ouverture automatique des .scn
    ;-- writeRegStr HKCR ".scn" "" "Sofa.scn" 
    ;-- writeRegStr HKCR "Sofa.scn\shell\open\command" "" '"$INSTDIR\bin\runSofa.exe %1"'
    ;-- writeRegStr HKCR "Sofa.scn\DefaultIcon" "" '"$INSTDIR\bin\runSofa.exe,0"'

    SetOutPath $INSTDIR\bin
    

  SectionEnd

  
  Section "sofaCUDA"  
  SectionIn 2
  SetOutPath $INSTDIR\bin
  File  "..\bin\sofaCUDA.exe"  
  
  CreateDirectory $INSTDIR\examples\CUDA
  SetOutPath $INSTDIR\examples\CUDA  
  File /r "..\examples\CUDA\*"  
  
  StrCpy $sofaCUDA_installed "yes"
  SectionEnd
  
SectionGroupEnd



;---------
SectionGroup "Tools"

  ;--------- autres executables
  Section "Modeler"
  SectionIn 1 2
    SetOutPath $INSTDIR\bin
    File "..\bin\modeler.exe"    
    
    CreateDirectory $INSTDIR\applications\projects\Modeler\preset
    SetOutPath $INSTDIR\applications\projects\Modeler\preset
    File /r "..\applications\projects\Modeler\preset\*"
    StrCpy $modeler_installed "yes"
  SectionEnd

  Section "generateRigid"
  SectionIn 2
    SetOutPath $INSTDIR\bin
    File "..\bin\GenerateRigid.exe"
  SectionEnd

  Section "meshConv"
  SectionIn 2
    SetOutPath $INSTDIR\bin
    File "..\bin\meshconv.exe"
  SectionEnd


SectionGroupEnd



;---------
SectionGroup "Documentation"

  ;-------- fichiers pdf
  Section "PDF manuals"
  SectionIn 1 2
    CreateDirectory $INSTDIR\doc
    SetOutPath $INSTDIR\doc
    File /r "..\doc\*.pdf"
    StrCpy $doc_installed "yes"
  SectionEnd

  ;--------- tous les exemples, et ce qu'il y a dans "share"
  Section "Main Demos"
  SectionIn 1 2
    CreateDirectory $INSTDIR\examples\Demos
    SetOutPath $INSTDIR\examples\Demos
    File /r "..\examples\Demos\*"
    StrCpy $examples_installed "yes" 
	SetOutPath $INSTDIR\share
	File /r "..\share\*"
  SectionEnd
  
  Section "Components examples"    
  SectionIn 1 2
	;---- Copy the Objects
	CreateDirectory $INSTDIR\examples\Objects
	SetOutPath $INSTDIR\examples\Objects
	File /r "..\examples\Objects\*"
	
	;---- Copy the .scn of the components
    CreateDirectory $INSTDIR\examples\Components
    SetOutPath $INSTDIR\examples\Components
    File /r "..\examples\Components\*"
    
    StrCpy $examples_installed "yes"
    StrCpy $examples_installed "yes" +3 0    
	SetOutPath $INSTDIR\share
	File /r "..\share\*"
  SectionEnd
  
  
  ;--Section "C++ Scenes"  
  ;--SectionIn 2
  ;--  CreateDirectory $INSTDIR\applications\tutorials
  ;--  SetOutPath $INSTDIR\applications\tutorials
  ;--  File /r "..\applications\tutorials\*"    
  ;--  StrCpy $examples_installed "yes"
  ;--  StrCpy $examples_installed "yes" +3 0  
  ;--  SetOutPath $INSTDIR\share
  ;--  File /r "..\share\*"
  ;--SectionEnd

  Section "Benchmark"  
  SectionIn 2
    CreateDirectory $INSTDIR\examples\Benchmark
    SetOutPath $INSTDIR\examples\Benchmark
    File /r "..\examples\Benchmark\*"   
  SectionEnd
  
SectionGroupEnd



;--------- Cr�ation des raccourcis et menus
Section "Start Menu & Desktop Shortcuts"
SectionIn 1 2
  CreateShortCut "$DESKTOP\runSofa.lnk" "$INSTDIR\bin\runSofa.exe" "" "$INSTDIR\bin\runSofa.exe" 0
  

  CreateDirectory "$SMPROGRAMS\Sofa${VERSION} "
  CreateShortCut "$SMPROGRAMS\Sofa${VERSION}\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\Sofa${VERSION}\runSofa.lnk" "$INSTDIR\bin\runSofa.exe" "" "$INSTDIR\bin\runSofa.exe" 0
  
  StrCmp $modeler_installed "yes" 0 +3
  CreateShortCut "$DESKTOP\Modeler.lnk" "$INSTDIR\bin\Modeler.exe" "" "$INSTDIR\bin\Modeler.exe" 0
  CreateShortCut "$SMPROGRAMS\Sofa${VERSION}\Modeler.lnk" "$INSTDIR\bin\Modeler.exe" "" "$INSTDIR\bin\Modeler.exe" 0


  StrCmp $sofaCUDA_installed "yes" 0 +3
  CreateShortCut "$DESKTOP\sofaCUDA.lnk" "$INSTDIR\bin\sofaCUDA.exe" "" "$INSTDIR\bin\sofaCUDA.exe" 0
  CreateShortCut "$SMPROGRAMS\Sofa${VERSION}\sofaCUDA.lnk" "$INSTDIR\bin\sofaCUDA.exe" "" "$INSTDIR\bin\sofaCUDA.exe" 0

  ;-- verifi� si la doc a �t� install�e. Met le raccourci si oui
  ;-- StrCmp $doc_installed "yes" 0 +2
  ;-- CreateShortCut "$SMPROGRAMS\Sofa${VERSION}\documentation directory.lnk" "$INSTDIR\doc" "" "$INSTDIR\doc" 0

SectionEnd




;--------------------------------

;---------- Uninstaller

Section "Uninstall"

  ;--- dabord on v�rifie les droit de l'utilisateur : si admin->install pour tout le monde si simple user->install que pour lui
  Pop $0
  UserInfo::GetAccountType
  Pop $1
  UserInfo::GetOriginalAccountType
  Pop $2
  StrCmp $1 "Admin" 0 +3
      SetShellVarContext all
    Goto done
  StrCmp $1 "Power" 0 +3
      SetShellVarContext all
    Goto done
  StrCmp $1 "User" 0 +3
      SetShellVarContext current
    Goto done
  StrCmp $1 "Guest" 0 +3
      SetShellVarContext current
    Goto done
  MessageBox MB_OK "Unknown error"
  Goto done

  done:
  
  ; suppression des cl�s de registre
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Sofa${VERSION}"
  DeleteRegKey HKLM SOFTWARE\Sofa

  ; liste des fichiers � d�sinstaller
  Delete $INSTDIR\bin\*.*
  Delete $INSTDIR\*.*

  SetOutPath $PROGRAMFILES
  RMDir /r $INSTDIR

  ; supression des raccourcis
  Delete "$SMPROGRAMS\Sofa${VERSION}\*.*"
  RMDir /r $SMPROGRAMS\Sofa${VERSION}
  Delete "$DESKTOP\runSofa.lnk"
  Delete "$DESKTOP\Modeler.lnk"
  Delete "$DESKTOP\sofaCUDA.lnk"


SectionEnd
