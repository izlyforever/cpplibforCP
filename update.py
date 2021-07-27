# python update.py
import os, sys
if len(sys.argv) < 2:
	os.system('git commit -a -m "update"')
else:
	os.system('git commit -a -m "' + ' '.join(sys.argv[1:]) + '"')
os.system('git push origin master')


# clear git history(uncomment select all copy):

# rm -rf .git 
# git init
# git add -A 
# git commit -m "reinit" 
# git remote add origin git@github.com:izlyforever/cpplibforCP.git
# git branch -M main
# git push -f -u origin main