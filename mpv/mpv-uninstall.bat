@echo off
cls
echo https://github.com/rossy/mpv-install
echo Files: mpv-install.bat, mpv.uninstall.bat, mpv-document.ico
echo.
echo Copyright 2015 James Ross-Gowan (rossymiles@gmail.com)
echo Permission to use, copy, modify, and/or distribute this software for any
echo purpose with or without fee is hereby granted, provided that the above
echo copyright notice and this permission notice appear in all copies.
echo.
echo THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
echo REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
echo FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
echo INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
echo LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
echo OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
echo PERFORMANCE OF THIS SOFTWARE.
echo.

setlocal enableextensions enabledelayedexpansion
path %SystemRoot%\System32;%SystemRoot%;%SystemRoot%\System32\Wbem

:: Unattended install flag. When set, the script will not require user input.
set unattended=no
if "%1"=="/u" set unattended=yes

:: Make sure the script is running as admin
call :ensure_admin

:: Delete "App Paths" entry
reg delete "HKLM\SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\mpv.exe" /f >nul

:: Delete HKCR subkeys
set classes_root_key=HKLM\SOFTWARE\Classes
reg delete "%classes_root_key%\Applications\mpv.exe" /f >nul
reg delete "%classes_root_key%\SystemFileAssociations\video\OpenWithList\mpv.exe" /f >nul
reg delete "%classes_root_key%\SystemFileAssociations\audio\OpenWithList\mpv.exe" /f >nul

:: Delete AutoPlay handlers
set autoplay_key=HKLM\SOFTWARE\Microsoft\Windows\CurrentVersion\Explorer\AutoplayHandlers
reg delete "%autoplay_key%\Handlers\MpvPlayDVDMovieOnArrival" /f >nul
reg delete "%autoplay_key%\EventHandlers\PlayDVDMovieOnArrival" /v "MpvPlayDVDMovieOnArrival" /f >nul
reg delete "%autoplay_key%\Handlers\MpvPlayBluRayOnArrival" /f >nul
reg delete "%autoplay_key%\EventHandlers\PlayBluRayOnArrival" /v "MpvPlayBluRayOnArrival" /f >nul

:: Delete "Default Programs" entry
reg delete "HKLM\SOFTWARE\RegisteredApplications" /v "mpv" /f >nul
reg delete "HKLM\SOFTWARE\Clients\Media\mpv\Capabilities" /f >nul

:: Delete all OpenWithProgIds referencing ProgIds that start with io.mpv.
for /f "usebackq eol= delims=" %%k in (`reg query "%classes_root_key%" /f "io.mpv.*" /s /v /c`) do (
	set "key=%%k"
	echo !key!| findstr /r /i "^HKEY_LOCAL_MACHINE\\SOFTWARE\\Classes\\\.[^\\][^\\]*\\OpenWithProgIds$" >nul
	if not errorlevel 1 (
		for /f "usebackq eol= tokens=1" %%v in (`reg query "!key!" /f "io.mpv.*" /v /c`) do (
			set "value=%%v"
			echo !value!| findstr /r /i "^io\.mpv\.[^\\][^\\]*$" >nul
			if not errorlevel 1 (
				echo Deleting !key!\!value!
				reg delete "!key!" /v "!value!" /f >nul
			)
		)
	)
)

:: Delete all ProgIds starting with io.mpv.
for /f "usebackq eol= delims=" %%k in (`reg query "%classes_root_key%" /f "io.mpv.*" /k /c`) do (
	set "key=%%k"
	echo !key!| findstr /r /i "^HKEY_LOCAL_MACHINE\\SOFTWARE\\Classes\\io\.mpv\.[^\\][^\\]*$" >nul
	if not errorlevel 1 (
		echo Deleting !key!
		reg delete "!key!" /f >nul
	)
)

echo Uninstalled successfully
if [%unattended%] == [yes] exit 0
pause
exit 0

:die
	if not [%1] == [] echo %~1
	if [%unattended%] == [yes] exit 1
	pause
	exit 1

:ensure_admin
	:: 'openfiles' is just a commmand that is present on all supported Windows
	:: versions, requires admin privileges and has no side effects, see:
	:: https://stackoverflow.com/questions/4051883/batch-script-how-to-check-for-admin-rights
	openfiles >nul 2>&1
	if errorlevel 1 (
		echo This batch script requires administrator privileges. Right-click on
		echo mpv-uninstall.bat and select "Run as administrator".
		call :die
	)
	goto :EOF
