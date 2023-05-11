from sys import *
from time import sleep
import multiprocessing as mp
mp.set_start_method('fork')

def play(connection, value):
   
    if value==0:
        print(f'First player serve')
    else:
        print(f'First player hit')
        
    value +=1
    
    
    connection.send(value)
 

def pingpong(connection, serve):
    
    if serve:
        play(connection, 0)
 
    while True:
        value = connection.recv()
        print(f'Second player hit \n')
        print(value)
        if value ==8:
            break
        
        play(connection, value)
    print("Game is over")
            

        
            
        
    
 

   
conn1, conn2 = mp.Pipe()
   
player1 = mp.Process(target=pingpong, args=(conn1,True))
player2 = mp.Process(target=pingpong, args=(conn2,False))
    
player1.start()
player2.start()
    
player1.join()
player2.join()