import tkinter as tk
from tkinter import *
import subprocess
from subprocess import Popen, PIPE
import webbrowser


def callback(): 
	n = en.get()
	t1 = et1.get()
	t2 = et2.get()
	l = el.get()
	k = ek.get()
	file = open("input.txt","w")
	file.write(n)
	file.write("\n")
	file.write(t1)
	file.write("\n")
	file.write(t2)
	file.write("\n")
	file.write(l)
	file.write("\n")
	file.write(k)
	file.close()
	subprocess.call(["g++", "./1Dheat.cpp"])
	tmp=subprocess.call("a.exe")
	webbrowser.open("output.txt")

r = tk.Tk()
r.title('1-D heat problem')
Label(r, text='Number of nodes').grid(row=0)
en = Entry(r)
Label(r, text='Temperature at first node').grid(row=1)
et1 = Entry(r)
Label(r, text='Temperature at last node').grid(row=2)
et2 =  Entry(r)
Label(r, text='Length of the rod').grid(row=3)
el = Entry(r)
Label(r, text='Value of k for the material').grid(row=4)
ek = Entry(r)

en.grid(row = 0, column = 1)
et1.grid(row = 1, column = 1)
et2.grid(row = 2, column = 1)
el.grid(row = 3, column = 1)
ek.grid(row = 4, column = 1)

button_submit = tk.Button(r, text='submit',bd = 4, bg = 'lightgreen', width=25, command=callback)
button_close = tk.Button(r, text='Close',bd = 4, bg = 'red', width=25, command=r.destroy)
button_submit.grid(row = 6, column = 1)
button_close.grid(row = 6, column = 2)
r.mainloop()
