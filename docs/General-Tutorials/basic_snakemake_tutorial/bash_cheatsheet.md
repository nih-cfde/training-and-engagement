# bash command cheatsheet

Keyboard shortcuts for terminal:

- use `tab` key to tab-complete command, file, and directory names

- use up/down keys to go to previous/current terminal commands

- `Ctrl z` to stop commands in progress

Bash | Description
--- | ---
`pwd` | path to current directory
`cd` | change directory
`ls` | list. `ls -lht` will list directory contents in human-readable (`-h`) file sizes and with time stamp (`-t`)
`history` | shows history of all commands you've entered. Note: the default command is `history 1` on zsh.
`less <file name>` | view file contents, scroll with up/down keys. exit view by typing `q`
\| | a pipe that connects commands
`history \| less` | view history
`nano -ET4 <file name>` | opens up the text editor `nano`, tab key = 4 spaces
`grep` | search files using regular expressions to pattern match. `grep rule Snakefile` will output all the lines in the file called Snakefile with the word 'rule' in them
`rm <file name>` | removes files forever, be careful!
`clear` | clear terminal window screen
`PS1='$ '` | change prompt symbol, e.g., to make it shorter. Note that if you are using a bash terminal, the prompt will end in an `$`. If you are using a zsh terminal, the prompt will end in an `%` by default.

Text editing with nano | description
--- | ---
`Ctrl k` | cut text
`Ctrl u` | paste text
`Ctrl _`, enter line number, return key| go to specific line number
`Ctrl e` | go to end of line
`Ctrl v` | scroll down page
`Ctrl o`, edit file name if needed, return key| saving changes (Write Out)
`Ctrl x` | exit nano

`#` in text files = comments; these lines of text are not interpreted as code
