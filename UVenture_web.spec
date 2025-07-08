# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['UVenture_web.py'],
    pathex=[],
    binaries=[],
    datas=[('static\\settings.txt', 'static'), ('static\\settings_default.txt', 'static'), ('static\\uventure_icon.png', 'static'), ('static\\help_page_contents.txt', 'static'), ('static\\taskstorage.txt', 'static'), ('templates\\*', 'templates'), ('UVenture_web.py', './'), ('UVenture.py', './'), ('UVenture_peakdetection.py', './'), ('MS_functions.py', './'), ('peakdetection_funcs.py', './')],
    hiddenimports=['scipy.stats._distn_infrastructure', 'scipy.special._ufuncs', 'scipy.stats', 'scipy.signal', 'scipy.stats.distributions'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['PyQt5'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='UVenture_web',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='UVenture_web',
)
