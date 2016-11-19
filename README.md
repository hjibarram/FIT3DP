# FIT3DP
First you need to define as executable file 

chmod a+x *.py dentro de la carpeta FIT3DP), después definir las variables de sistema 
export FIT3DP_PATH="directorio de la carpeta FIT3DP"
export PATH="$FIT3DP_PATH:$PATH"
después de esto en teoría deberia de correr en cualquier dierctorio, por ejemplo para analizar un cubo que el usuario desee usando el comando
ana_single.py manga-xxxx-yyyy
Nota, hay que definir las variables DIR_DATA,DIR_DATA_OUT,DIR_PLOTS, del archivo "ana_single.py" acorde de donde se encuentran los archivos a analizar y donde se quieren guardar