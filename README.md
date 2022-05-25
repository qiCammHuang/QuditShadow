# QuditShadow
A MATLAB+Python demo for classical shadow on Qudits.






## Python related usage

Before running ```main_highdim.m```, make sure that your MATLAB has successfully loaded a python executable. You could check by using the ```pyversion``` command. For example,

```matlab
>> pyversion

       version: '3.7'
    executable: '/Users/user/opt/anaconda3/envs/classicalShadow/bin/python3.7'
       library: '/Users/user/opt/anaconda3/envs/classicalShadow/lib/libpython3.7m.dylib'
          home: '/Users/user/opt/anaconda3/envs/classicalShadow'
      isloaded: 1
```


If changes are made to python files, MATLAB would not immediately execute the new programme. Instead, in our case, do

```matlab
% Reload Python module
clear classes
obj = py.importlib.import_module('derandTools');
py.importlib.reload(obj);
```

so that MATLAB can run the renewed python files.


