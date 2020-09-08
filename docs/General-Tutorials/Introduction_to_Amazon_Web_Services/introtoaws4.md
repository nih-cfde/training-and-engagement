---
layout: page
title: Connect to an Instance
---
## Step 3: Connect to your Instance

=== "Windows Instruction"

    Ok, so you've created a running computer. How do you get to it?

    The main thing you'll need is the network name of your new computer. To retrieve this, go to the instance view and click on the instance, and find the "Public DNS." This is the public name of your computer on the Internet.

    Copy this name, and connect to that computer with ssh under the username 'ubuntu,' as follows.

    ## Install MobaXTerm

    MobaXTerm is a terminal for Windows with an X11 server, a tabbed SSH client and several other network tools for remote computing (VNC, RDP, telnet,rlogin)
    MobaXterm brings all the essential Unix commands to Windows desktop, in a single portable exe file which works out of the box.

    First, download mobaxterm and run it.

    Download Link: https://download.mobatek.net/2032020060430358/MobaXterm_Portable_v20.3.zip

    ## Start a new session

    ![mobaxterm1](../../images/mobaxterm_1.PNG)

    ![Remote Host](../../images/mobaxterm_2.PNG)

     ## Fill in session settings

    ## Specify the session key

    Put in your hostname (should be ec2-XXX-YYY-AAA.compute-1.amazon.aws.com) and enter 'ubuntu'

    ![Hostname](../../images/mobaxterm_3.PNG)

    Copy the amazon.pem file into your primary hard disk (generally C:) and then put the full path into it.

    ![Private Key](../../images/mobaxterm_3_2.png)

    ## Click OK

    If you see this screen and the line ubuntu@ip-###-##-#-##:~$,
    this means your instance computer is ready.

    ![Ubuntu Terminal](../../images/mobaxterm_4.PNG)

    ## Transferring Files Using Mobaxterm
    - With Mobaxterm, you can transfer files between your local computer and your remote instance by simply dragging and dropping files between MoBaxterm's SCP tab (located on the left-hand side of the MobaXTerm window) and your local computer's file explorer.
    ![SCP Tab](../../images/Mobaxterm_transfer1.PNG)
    
    ![Transfer File](../../images/Mobaxterm_transfer2.PNG)

        


    ## Step 5: Terminating the Instance

    Once you have completed your tasks and are sure you do not need the instance any longer, you may terminate the instance by returning to AWS Management Console

    !!! Warning
        
        Terminating an instance will erase all the work you have saved on the instance! Be sure to save your work on your local computer or on a cloud storage drive.

    - Click on 'Services'
    - Click 'EC2' 
    - Click 'Instance' on the left hand side bar menu and it should bring you to the list of running instances on your account.
    - Click on the instance you would like to terminate
    - Click 'Actions'
    - Click 'Instance State'
    - Select 'Terminate'

    !!! Note
        If you just close the terminal, it doesn't mean the instance has been terminated. It will keep going until you terminate at the AWS console.



        ![Terminate](../../images/Terminate.png)



=== "Mac OS"
    
    OK, so you've created a running computer. How do you get to it?

    You will need the network name of your new computer. TO retreive this, go the instance view and click on the instance, and find the "Public DNS." This is the public name of your computer on the Internet.

    Copy this name, and connect to that computer with ssh under the username 'ubuntu', as follows

    First, find your private key file; it is the .pem file you downloaded when starting up your EC2 instance. It should be in your Downloads folder. Move it onto your desktop.

    Next start Terminal and type:
    ```
    chmod og-rwx ~/Desktop/######.pem
    ```
    where ###### is the name of the .pem file we downloaded earlier

    Then type;
    ```
    ssh -i ~/Desktop/######.pem ubuntu@ ec2-???-???-???-???.compute-1.amazonaws.com
    ```
    where ec2-???-???-???-???.compute-1.amazonaws.com is the Public DNS we copied earlier.

    ## Step 4: Transferring Files

    ## Copying files EC2 to Local Computer

    - To use `scp` (secure copy) with a key pair use the following command:
    
    ```
    scp -i /directory/to/amazon.pem ubuntu@ec2-xx-xxx-xxx.compute-1.amazonaws.com:path/to/file /your/local/directory/files/to/download
    ```

    - You may also download a file from EC2, download folder by archiving it

    zip -r squash.zip /your/ec2/directory/

    You can download all archived files from ec2 by typing the following command
    ```
    scp - i/directory/to/amazon.pem ubuntu@ec2-xx-xx-xxx-xxx.compute-1.amazonaws.com:~/* /your/local/directory/files/to/download
    ```

    ## Copying files from local to EC2

    - Your private key must not be publicly visible. Run the following command so that only the root user can read the file.
    ```
     chmod 400 amazon.pem
    ```
    To use `scp` with a key pair use the following command:
    ```
    scp -i /directory/to/amazon.pem /your/local/file/to/copy ubuntu@ec2-xx-xx-xxx-xxx.compute-1.amazonaws.com:path/to/file

    !!! Note:

        You need to make sure that the user "user" has the permission to write in the target directory. In this example, if ~/path/to/file was created by user "user", it should fine.



    ## Step 5: Terminating the Instance

    Once you have completed your tasks and are sure you do not need the instance any longer, you may terminate the instance by returning to AWS Management Console

    !!! Warning
        
        Terminating an instance will erase all the work you have saved on the instance! Be sure you save all of your work that you would like to keep on your local computer or cloud storage drives.

    - Click on 'Services'
    - Click 'EC2' 
    - Click 'Instance' on the left hand side bar menu and it should bring you to the list of running instances on your account.
    - Click on the instance you would like to terminate
    - Click 'Actions'
    - Click 'Instance State'
    - Select 'Terminate'

    !!! Note
        If you just close the terminal, it doesn't mean the instance has been terminated. It will keep going until you terminate at the AWS console.



        ![Terminate](../../images/Terminate.png)

