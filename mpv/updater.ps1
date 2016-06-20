$url = "https://bitbucket.org/rorgoroth/mpv-for-windows/downloads/mpv-latest-Win64.zip"
$output = "mpv-latest-Win64.zip"

Invoke-WebRequest -Uri $url -OutFile $output
$shell_app= new-object -com shell.application
$filename = "mpv-latest-Win64.zip"
$zip_file = $shell_app.namespace((Get-Location).Path + "\$filename")
$destination = $shell_app.namespace((Get-Location).Path)
$destination.Copyhere($zip_file.items())

Remove-Item -Force .\mpv-latest-Win64.zip
