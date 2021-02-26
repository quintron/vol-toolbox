cd build
cmake --build . --config Release --target install
ctest -C Release
cd..

pip install -r requirements.txt
python setup.py develop

pip install pytest
pytest tests