# python update.py
import os, sys
os.system('git add .') 
if len(sys.argv) < 2:
	os.system('git commit -m "update"')
else:
	os.system('git commit -m "' + ' '.join(sys.argv[1:]) + '"')
os.system('git push origin main')


# clear git history:

# sudo rm -rf .git 
# git init
# git add -A 
# git commit -m "init" 
# git remote add origin git@github.com:izlyforever/cpplibforCP.git
# git branch -M main
# git push -f -u origin main