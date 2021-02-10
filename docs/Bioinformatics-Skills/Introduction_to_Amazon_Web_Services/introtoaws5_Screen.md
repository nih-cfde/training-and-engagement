# The Screen Command

Screen or GNU Screen is a terminal multiplexer. You can start a terminal session and then open multiple screens inside that session. Processes running in Screen will continue to run when their window is not visible, even if you get disconnected. The goal of this tutorial is for users to learn how to run and switch between multiple screen sessions

=== "Learning Objectives"
    - Learn to run multiple screen sessions using Linux Screen command
    - Learn to switch between multiple screen sessions
    - Learn to recover a running session in case a sessions disconnects
=== "Est. Time"
    - 20 minutes
=== "Prerequisites"
    - Run an AWS Instance: Ubuntu Pro 20.04 LTS
    - Connect to an AWS instance. Windows users will require MobaXterm (for windows users)


## Installing the Screen

To install screen, run the following command:

=== "AWS Instance Code"
    ```
    sudo apt update
    sudo apt-get install screen
    ```

=== "Expected Output"
    ```
    ubuntu@ip-172-31-7-6:~$ sudo apt-get install screen
    Reading package lists... Done
    Building dependency tree       
    Reading state information... Done
    screen is already the newest version (4.8.0-1).
    0 upgraded, 0 newly installed, 0 to remove and 7 not upgraded.
    ```

The first line of code updates your instance to all latest software configurations. The second line of code does the actual installing.

## Running Screen

To help yourself tell the different terminal screen apart, type this command on the original terminal window (before you try out screen):

=== "AWS Instance Code"

    ```
    echo "this is the original terminal window"
    ```


Then type the command `screen` into the same AWS instance to start a new screen session.

=== "AWS Instance Code"
    ```
    screen
    ```
=== "Expected Output"
    ```
    GNU Screen version 4.08.00 (GNU) 05-Feb-20

    Copyright (c) 2018-2020 Alexander Naumov, Amadeusz Slawinski
    Copyright (c) 2015-2017 Juergen Weigert, Alexander Naumov, Amadeusz Slawinski
    Copyright (c) 2010-2014 Juergen Weigert, Sadrul Habib Chowdhury
    Copyright (c) 2008-2009 Juergen Weigert, Michael Schroeder, Micah Cowan,
    Sadrul Habib Chowdhury
    Copyright (c) 1993-2007 Juergen Weigert, Michael Schroeder
    Copyright (c) 1987 Oliver Laumann

    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free Software
    Foundation; either version 3, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with
    this program (see the file COPYING); if not, see

                  [Press Space for next page; Return to end.]
    ```
Press Space (twice) or Return to get to the command prompt. You are now on a new screen

## Using screen

Let's run a program in the new screen window to test it out.

=== "AWS Instance Code"
    ```
    top
    ```
The `top` command is used to show the Linux processes. It provides a dynamic real-time view of the running system. Usually, this command shows the summary information of the system and the list of processes or threads which are currently managed by the Linux Kernel.

While the top command is still running, create another screen by clicking ++ctrl+a+c++. You should see a new blank terminal.

Let's run yet another command on the blank terminal:

=== "AWS Instance Code"
    ```
    ping google.com
    ```
=== "Expected Output"
    ```
    PING google.com (172.217.4.206) 56(84) bytes of data.
    64 bytes from ord37s19-in-f14.1e100.net (172.217.4.206): icmp_seq=1 ttl=100 time=17.4 ms
    64 bytes from ord37s19-in-f14.1e100.net (172.217.4.206): icmp_seq=2 ttl=100 time=17.4 ms
    64 bytes from ord37s19-in-f14.1e100.net (172.217.4.206): icmp_seq=3 ttl=100 time=17.5 ms
    64 bytes from ord37s19-in-f14.1e100.net (172.217.4.206): icmp_seq=4 ttl=100 time=17.5 ms
    .....
    ```

A ping test is a method of checking if the computer is connected to a network. It is used to ensure that a host computer which your computer tries to access is operating. It is run for troubleshooting to know connectivity as well as response time.

Please note that you wow have three running screens:
 Screen 1) The original ssh terminal you saw at login. You typed the echo command in it.
 Screen 2) A screen that's running the `top` command.
 Screen 3) A screen that's running the `ping google.com` command.

### Switching between screens

 To switch between the two new screens, i.e. screen 2 and screen 3, type ++ctrl+a+p++. You cannot toggle to the original terminal screen (i.e. screen 1) with this shortcut.

### Detaching screens

 To detach a screen session and return to your original SSH terminal, type ++ctrl+a+d++. You will be taken back to the terminal window with the `echo` command:

![](./images-aws/original_ssh_terminal.png)

To list your current screen sessions type:

=== "AWS Instance Code"
    ```
    screen -ls
    ```
=== "Expected Output"
    ```
    ubuntu@ip-172-31-7-6:~$ screen -ls
    There is a screen on:
	   2683.pts-0.ip-172-31-7-6	(02/09/21 19:41:19)	(Detached)
     1 Socket in /run/screen/S-ubuntu.
    ```

You can use the screen id to reconnect to your screen. Like this:

=== "AWS Instance Code"
    ```
    screen -r 2683.pts-0.ip-172-31-7-6
    ```
You should see Screen 2 that you previously created. Once again you can toggle between screen 2 and screen 3 by typing ++ctrl+a+p++


## Quitting screens

You can quit the `top` output by typing `q` or ++ctrl+c++

To end a screen session, toggle into the session and type:

=== "AWS Console Code"
    ```
    exit
    ```

If no other screen sessions are open, you will fall back to the original SSH terminal (screen 1). If another screen session is open, you will enter that screen session. You can type `exit` until all screen sessions are closed.

!!! Warning
    Typing `exit` into a screen permanently closes that screen session and all the analyses that are conducted in it



## Screen Cheat Sheet

 Command | Description
--------|---------
++ctrl+a+c++ | Creates a new screen session so that you can use more than one screen session at once
++ctrl+a+n++ | Switches to the next screen session if you have more than one running screen
++ctrl+a+p++ | Switches to the previous screen session if you have more than one running screen
++ctrl+a+d++ | Detaches a screen session without killing the processes running in it
`exit`| Kills a screen session permanently

#### Getting in
Command | Description
--------|---------
`screen -S <name>` | Start a new screen session with session name
`screen -ls` | List running sessions/screens
`screen -x` | Attach to a running session
`screen -r <name>` |Attach to a running session with a name

#### Getting Out
Command | Description
--------|---------
`screen -d <name>` | Detach a running session
++ctrl+a+d++ | Detaches a screen session without killing the processes running in it
++ctrl+a++ then ++ctrl+d++| Detach and logout (quick exit)

#### Help
Command | Description
--------|---------
++ctrl+a++ | See help
++ctrl+a+c++ | Create new window
++ctrl+a+n++ or ++ctrl+a++ <space> | Change to next window in list
++ctrl+a+p++ or ++ctrl+a++ <backspace> | Change to previous window in list

video tutorial link
https://video.ucdavis.edu/media/AWS_Ubuntu_Screen_v2/1_lkmj4ful/166161802
