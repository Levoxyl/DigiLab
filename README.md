# 🛠 Requirments
# Linux (WSL)
## Externally managed enviroment

```

cd ~/DigiLab

python3 -m venv venv

source venv/bin/activate

``` 

## Dependencies

### BioPython & dotenv  

~~~

pip install python-dotenv biopython

~~~

*Note: Select python interpreter as venv so you don't get dependency errors

# Windows


## 1. Initialize Python Virtual Environment
```
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install PyQt6 PyQt6-WebEngine biopython python-dotenv
```

## 2. Setup Frontend Dependencies
```
$env:Path += ";<YOUR_NODE_PATH>" # e.g., D:\Code\Tools\Node
cd App\front\themes\modern
npm install
```

# Powershell Running  

## Run frontend
```
cd App\front\themes\modern
npm run dev
```

## Run backend/UI launcher
```
.\venv\Scripts\Activate.ps1
python App\front\launcher.py
```

