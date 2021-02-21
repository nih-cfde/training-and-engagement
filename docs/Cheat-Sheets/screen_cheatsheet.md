# Screen Cheat Sheet

 Command | Description
--------|---------
++ctrl+a+c++ | Creates a new screen session so that you can use more than one screen session at once
++ctrl+a+n++ | Switches to the next screen session if you have more than one running screen
++ctrl+a+p++ | Switches to the previous screen session if you have more than one running screen
++ctrl+a+d++ or ++ctrl+a++ then ++ctrl+d++ | Detaches a screen session without killing the processes running in it
++ctrl+a++ then ++ctrl+a++ | Switches between screens
`exit`| Kills a screen session permanently

## Getting in
Command | Description
--------|---------
`screen -S <name>` | Start a new screen session with session name
`screen -ls` | List running sessions/screens
`screen -r` | Attach to a running session
`screen -r <name>` |Attach to a running session with a name

## Getting Out
Command | Description
--------|---------
`screen -d <name>` | Detach a running session
++ctrl+a+d++ | Detaches a screen session without killing the processes running in it
++ctrl+a++ then ++ctrl+d++| Detach and logout (quick exit)
screen -S <screen#/name> -X quit | Delete a screen while in detached state

## Toggling
Command | Description
--------|---------
++ctrl+a+c++ | Create new window
++ctrl+a+n++ or ++ctrl+a++ <space> | Change to next window in list
++ctrl+a+p++ or ++ctrl+a++ <backspace> | Change to previous window in list

## Help
Command | Description
--------|---------
screen -h | See help
