@echo off
cd /d "%~dp0"
echo Installing FastAPI dependencies...
C:\Users\webst\AppData\Local\Programs\Python\Python314\python.exe -m pip install fastapi uvicorn python-multipart --quiet
echo.
echo Starting Chemical Engineering Solver Web API...
echo.
echo Web Interface: http://localhost:8000
echo API Docs: http://localhost:8000/docs
echo.
C:\Users\webst\AppData\Local\Programs\Python\Python314\python.exe api.py
pause
