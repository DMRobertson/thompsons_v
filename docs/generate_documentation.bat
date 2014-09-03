@Echo Off

echo.Running apidoc.
sphinx-apidoc -o . -e -F -H "Thompson's V" -A "David Robertson" ..

echo.Cleaning up build directory.
call make clean

echo.Running tests.
call make doctest

echo.Building html.
call make html
