cd docs/
sphinx-apidoc -force -o ./api/fenicsx/. ../models/fenicsx/

rm api/fenicsx/modules.rst
cp -r ../Tutorials . # needed since sphinx doesn't include parent directoris, the copied directory is not pushed

make clean
make html
