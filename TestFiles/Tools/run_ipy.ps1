# Draait een IronPython-script met de IronPython 2.7.9 van Mecway, zonder Mecway
# (controle van Python 2.7-compatibiliteit en offscreen GUI-tests).
# gebruik: powershell -ExecutionPolicy Bypass -File run_ipy.ps1 -Script <pad.py> -ScriptArgs a,b
param([string]$Script, [string[]]$ScriptArgs)
$mw = "C:\Program Files\Mecway"
foreach ($d in "Microsoft.Scripting","Microsoft.Dynamic","IronPython","IronPython.Modules") {
    Add-Type -Path (Join-Path $mw "$d.dll")
}
$engine = [IronPython.Hosting.Python]::CreateEngine()
$paths = $engine.GetSearchPaths()
$paths.Add("$mw\IronPython.2.7.9\Lib")
$paths.Add((Split-Path $Script))
$engine.SetSearchPaths($paths)
$out = New-Object System.IO.MemoryStream
$engine.Runtime.IO.SetOutput($out, [System.Text.Encoding]::UTF8)
$engine.Runtime.IO.SetErrorOutput($out, [System.Text.Encoding]::UTF8)
$sys = [IronPython.Hosting.Python]::GetSysModule($engine)
$argv = New-Object IronPython.Runtime.List
$argv.append($Script)
foreach ($a in $ScriptArgs) { $argv.append($a) }
$sys.SetVariable("argv", $argv)
try {
    $engine.ExecuteFile($Script) | Out-Null
} catch {
    $eo = $engine.GetService([Microsoft.Scripting.Hosting.ExceptionOperations])
    Write-Output ("FOUT: " + $eo.FormatException($_.Exception.InnerException))
}
[System.Text.Encoding]::UTF8.GetString($out.ToArray())
