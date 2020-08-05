---
layout: page
title: Enabling User Choices with Javascript
---

<script language="javascript" type="text/javascript">
function set_page_view_defaults() {
    document.getElementById('div_markdown').style.display = 'block';
    document.getElementById('div_html').style.display = 'none';
};

function change_content_by_syntax(form_control){
    if (!form_control || document.getElementById(form_control).value == 'value_markdown') {
        set_page_view_defaults();
    } else if (document.getElementById(form_control).value == 'value_html') {
        document.getElementById('div_markdown').style.display = 'none';
        document.getElementById('div_html').style.display = 'block';
    } else {
        alert("Error: Missing platform value for 'change_content_by_syntax()' script!");
    }
}

window.onload = set_page_view_defaults;
</script>

# Adding user choosable content using Javascript

There are many instances in a tutorial where you might want to give the user an option
to only see parts of the content. For instance, you might have instructions for how to
perform a task by operating system, and therefore would like to show the user only the
content relevant to their OS. This is useful, because many errors in workshops are due to
users being confused about which instructions apply to them, accidentally following the
instructions for all OS's instead of just one, and other similar issues.

One way to limit what a user sees is by running scripts in the page that give the user the
ability to dynamically change the content. In this tutorial, I'll show you how to use
Javascript to make a simple content chooser that changes the display based on the users
choice of OS.

You will need:

  - a raw web page in markdown or HTML with sections you want to be user choosable
  - a plain text editor
  - about 30 minutes


This tutorial can be used for websites written in markdown or pure HTML.

**Please select the syntax you will be using: <select id="id_syntax" name="platformlist" onchange="change_content_by_syntax('id_syntax');return false;"><option value="value_markdown" id="id_markdown" selected> markdown </option><option value="value_html" id="id_html" > HTML </option>></select>**


# Overview of the process

There are three steps to creating the javascript code in a prebuilt page:

1. Each chunk of page that should be hidden/shown by the user choice is wrapped in a bit of
HTML. These are called 'blocks'.
2. Add the script code logic to the raw page file
3. At the point of the page where you want your user to have the choice, insert the
code for the dropdown box


## Create blocks

In your current document, identify which parts of the document you want
to always be seen, and which you want to respond to the users choice. For each
section that should respond to the user. Add a `<div>` element around the section.
A `<div>` element is a bit of HTML that works as a divider, and we're using it
to section off this bit of page so our javascript will know what the boundaries
of the sections are. The most basic one would look like this:

```
<div>

My Windows only text.

</div>
```

As with all HTML elements, there is a beginning and an end. `<div>` means 'start'
the divider, `</div>` means 'end the divider'.

We've divided this section out, but in order for our code to identify it, it will
need a name. So we want to add some extra information into our `<div>`. Since my
sentence tells me that this part of my page is about Windows, let's make this
`<div>`s ID "div_windows":

```
<div id="div_windows">

My Windows only text.

</div>
```

We can also give it a style, that is, how do we want it to look or behave. For
instance, if we were inserting a table, we'd want the table style. Since this
section is just text, we'll give it a generic style:

```
<div id="div_windows" style="display:block">

My Windows only text.

</div>
```

<div id="div_markdown" style="display:block" markdown="1">

You may have noticed that we're putting a lot of HTML code into a markdown
document. That's ok! Most markdown renderers just pass html through as they render,
so HTML will show up in your final page without any changes. However, that means
that your markdown renderer ignores *everything* inside HTML tags and simply passes
them through. Since we've wrapped our whole section in a `<div>`, that means we
have now lost all of our markdown formatting in this section. To get it back,
we need to tell the markdown renderer to read through our `<div>` and render
anything that looks like markdown:

```
<div id="div_windows" style="display:block" markdown="1">

An windows option of text

</div>
```

</div>

Continue adding `div` tags around each of your sections. Don't put anything
around sections that you want displayed at all times. Be sure to give different
sections different IDs that match their contents. You can make the IDs anything
you like, however as you're learning, I suggest adding 'div_' to the front of
each ID so that you remember what that variable is for.

When you're done, your file should look something like this:

<div id="div_markdown" style="display:block" markdown="1">

```

Text that is always displayed

<div id="div_windows" style="display:block" markdown="1">

An windows option of text

</div>

More text that is always displayed

<div id="div_unix" style="display:block" markdown="1">

A unix option of text

</div>

<div id="div_mac" style="display:block" markdown="1">

A mac option of text

</div>

Even more text that is always displayed

```
</div>

<div id="div_html" style="display:block" markdown="1">

```

Text that is always displayed

<div id="div_windows" style="display:block">

An windows option of text

</div>

More text that is always displayed

<div id="div_unix" style="display:block">

A unix option of text

</div>

<div id="div_mac" style="display:block">

A mac option of text

</div>

Even more text that is always displayed

```
</div>


## Add the script code to the raw page file

In our example text, we have three different blocks separated by `div`s:
"windows", "unix", and "mac".

That means that our logic script will have to differentiate between those three
blocks of the page. We'll use some simple if/else logic to tell the page which
parts it should show the user, and which it shouldn't, depending on what option
the user picks in a dropdown box. (We'll code the dropdown box in the next section).


<div id="div_markdown" style="display:block" markdown="1">

Usually, scripts go at the top of the page, just below any YAML header.

</div>

<div id="div_html" style="display:block" markdown="1">

This script code should end up in the 'header' block of your HTML.

</div>

Since we are writing all of this in HTML, it should not be surprising that we
need to wrap our script in tags. Since we want this to run as a script, we use
the script tag, and we know we're going to write our code in javascript so we
specify that as well:

```
<script language="javascript" type="text/javascript">

</script>

```

First, let's write the code for what the page should do by default, that is,
what should the page show before our user has made a selection?

To do that, we're going to write a small function. Let's make an empty function
called 'set_page_view_defaults':

```
<script language="javascript" type="text/javascript">

function set_page_view_defaults() {
}

</script>

```

!!! Note

    This function is not a default part of javascript, it, and the other functions in this
    tutorial are custom functions, which is why they start with 'function': it's a
    declaration of what the following code is. So, you can't look up these functions, they
    don't have help pages, they're exactly like any shell script or R script or Python
    script you've written, just in a language your less familiar with.

We want the function to either show or not show each of our sections, since this
is the default, it isn't even going to check to see what the dropdown says, it
will just display this view.

Let's say we decide we want the Windows section to be the default display.
That means we want 'div_windows' to have the style 'block'
as we previously specified, but we want the others to not display. To easily not
display them, let's change their style to 'none'. We call each one by asking for
the page element of our document, using it's ID:

```
<script language="javascript" type="text/javascript">

function set_page_view_defaults() {
    document.getElementById('div_win').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
}

</script>
```

That takes care of the default. Now we need a function that can change the page
based on user input. To start, we'll again make an empty function, lets call it
change_content_by_platform:

```
<script language="javascript" type="text/javascript">

function set_page_view_defaults() {
    document.getElementById('div_win').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
}

function change_content_by_platform(){
}

</script>
```

In our first function, there was no input, it is just what happens by default in
the page. To have a responsive function, it will have to receive input from
somewhere. Let's call that input 'form_control':

```
<script language="javascript" type="text/javascript">

function set_page_view_defaults() {
    document.getElementById('div_win').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
}

function change_content_by_platform(form_control){
}

</script>
```

Now, we can add the logic into our function. We'll use the same kind of logic we
used for setting the default. In fact, if the user sets the dropdown to "Windows"
we want the page to literally respond with the default view, so let's add that as
an 'if statement':

```
<script language="javascript" type="text/javascript">

function set_page_view_defaults() {
    document.getElementById('div_win').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
}

function change_content_by_platform(form_control){
    if (!form_control || document.getElementById(form_control).value == 'value_win') {
        set_page_view_defaults()
    }
}

</script>
```

Our logic now says "if there is no value in form_control then set the value of
form_control to 'value_win', and run the function that shows the default view."


What is 'value_win'? Well, we haven't written our dropdown code yet, but we know
that it will output a variable that we've decided to call 'form_control'. We also
know that we will want to program 'form_control' to have 3 possible values, one
for each OS. We could call those 3 values whatever we want, but while you're
learning, it's easiest to keep the names simple. So I'm calling mine 'value_win',
'value_unix', and 'value_mac' to go with my 3 'div_'s.


Now let's account for user input. We've already written an if statement for the
default that just calls our set_page_view_defaults function. We can write the
unix option in the same basic way: an if statement that checks the value of
form_control, and a list of what elements should be shown for a given value. Note
that since we already have an if, we're now using 'else if':

```
<script language="javascript" type="text/javascript">

function set_page_view_defaults() {
    document.getElementById('div_win').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
}

function change_content_by_platform(form_control){
    if (!form_control || document.getElementById(form_control).value == 'value_win') {
        set_page_view_defaults();
    } else if (document.getElementById(form_control).value == 'value_unix') {
        document.getElementById('div_windows').style.display = 'none';
        document.getElementById('div_unix').style.display = 'block';
        document.getElementById('div_mac').style.display = 'none';
}

</script>
```

=== "Exercise"

    Try adding another 'else if' statement that gives our function the logic
    for how to display the mac parts of the page.

=== "Solution"

    ```
    else if (document.getElementById(form_control).value == 'value_mac') {
         document.getElementById('div_windows').style.display = 'none';
         document.getElementById('div_unix').style.display = 'none';
         document.getElementById('div_mac').style.display = 'block';
    }
    ```

If we're going to put our script in a production webpage, its best to have at
least a bit of error handling. Otherwise, the page loads too slowly, it might
have a hard to diagnose problem where the logic script runs before the page loads
the dropdown code, so there is no 'form_control' variable to check. So, let's
tell our script what to do if there is no input variable at all:

```
<script language="javascript" type="text/javascript">
function set_page_view_defaults() {
    document.getElementById('div_windows').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
};

function change_content_by_platform(form_control){
    if (!form_control || document.getElementById(form_control).value == 'value_win') {
        set_page_view_defaults();
    } else if (document.getElementById(form_control).value == 'value_unix') {
        document.getElementById('div_windows').style.display = 'none';
        document.getElementById('div_unix').style.display = 'block';
        document.getElementById('div_mac').style.display = 'none';
   } else if (document.getElementById(form_control).value == 'value_mac') {
        document.getElementById('div_windows').style.display = 'none';
        document.getElementById('div_unix').style.display = 'none';
        document.getElementById('div_mac').style.display = 'block';
    } else {
        alert("Error: Missing platform value for 'change_content_by_platform()' script!");
    }
}

</script>

```

Now, instead of breaking, our script can print an alert message and close gracefully.

There is one more thing to add before our script is complete, and that is an
event handler. This is a little bit of code that makes your script responsive to
input. What we have already written makes our script run whenever a user makes
a choice, but the default is *no user input*. That means that our script currently
won't run until the user does something. What we actually want, is for it to run
every time the user changes the choice *and* to run once when the page is loaded,
so it can set the default:

```
<script language="javascript" type="text/javascript">
function set_page_view_defaults() {
    document.getElementById('div_windows').style.display = 'block';
    document.getElementById('div_unix').style.display = 'none';
    document.getElementById('div_mac').style.display = 'none';
};

function change_content_by_platform(form_control){
    if (!form_control || document.getElementById(form_control).value == 'value_win') {
        set_page_view_defaults();
    } else if (document.getElementById(form_control).value == 'value_unix') {
        document.getElementById('div_windows').style.display = 'none';
        document.getElementById('div_unix').style.display = 'block';
        document.getElementById('div_mac').style.display = 'none';
   } else if (document.getElementById(form_control).value == 'value_mac') {
        document.getElementById('div_windows').style.display = 'none';
        document.getElementById('div_unix').style.display = 'none';
        document.getElementById('div_mac').style.display = 'block';
    } else {
        alert("Error: Missing platform value for 'change_content_by_platform()' script!");
    }
}

window.onload = set_page_view_defaults;

</script>

```

!!! Tip

    If you will be using these logic functions frequently throughout your website, with the same
    options, you can add this part of the script to your HTML templates in the header and
    have it run automatically on all pages. Running scripts will slow down your website,
    so there is a tradeoff between this ease in coding and the speed of your user experience.



## Add the dropdown box

Now we need to find a logical location for your dropdown box. This should be
before the first instance of a `div`. Depending on the flow of your
document, it might make the most sense to put it right at the top, as the first
thing your user sees, or it might make more sense to put it after some explanation
as I did in this page.

Let's start with some text explaining what your user should do.

<div id="div_html" style="display:block" markdown="1">

You can add any HTML formatting tags, such as bolding here as well.

```
<b>Please select the platform you wish to use for this exercise:</b>
```

</div>

<div id="div_markdown" style="display:block" markdown="1">

This will be rendered with markdown, so you can use markdown formatting such
as bolding with asterisks.

```
**Please select the platform you wish to use for this exercise:**
```

</div>


To make the dropdown, we'll need to set up some variables in the HTML that can
send information to our function. This variable will be populated by the users
selection, so we use the `<select>` tag. As we did with our previous HTML tags,
we give this one an ID, and I have chosen a descriptive one to keep things simple:

First, let's give our variable a name and ID

<div id="div_html" style="display:block" markdown="1">

```
<b>Please select the platform you wish to use for this exercise:</b> <select id="id_platform" >

```
</div>

<div id="div_markdown" style="display:block" markdown="1">

```
**Please select the platform you wish to use for this exercise:**
<select id="id_platform" >
```
</div>

Since this is a dynamic piece of our page, we need to give the `<select>` a
directive as well. When it is changed, we want it to send that change to our
logic function:

<div id="div_html" style="display:block" markdown="1">

```
<b>Please select the platform you wish to use for this exercise:</b>

<select id="id_platform" onchange="change_content_by_platform('id_platform');return false;">
```
</div>

<div id="div_markdown" style="display:block" markdown="1">

```
**Please select the platform you wish to use for this exercise:**

<select id="id_platform" onchange="change_content_by_platform('id_platform');return false;">
```
</div>

So, this code will accept the users input, save that input value as the variable
'id_platform' and send that variable to be the input into our change_content_by_platform
function, at which point 'form_control' will take on that value. The 'return false'
here tells the `<select>` that once it has passed on the user input, it shouldn't
do anything else.

Finally, we need to tell the `<select>` what the actual dropdown menu should be.
We do this by telling it all of the options it should have using the `<option>`
tag *inside* of our `<select>`. Our option tag value *must* match the values we
put in our logic function. The ID can be anything, but again I'm using matching
values that start with 'id_'.This is the `<option>` code for the Mac choice:

<div id="div_html" style="display:block" markdown="1">

```
<b>Please select the platform you wish to use for this exercise:<b>

<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
</select>

```

</div>

<div id="div_markdown" style="display:block" markdown="1">

```
**Please select the platform you wish to use for this exercise:**

<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
</select>

```

</div>

Note that there is a bit of text between the opening of `<option>` and the end
`</option>`. The text inside the option block is the text that will appear to
the user in the dropdown box.

=== "Exercise"

    Write the `<option>` values for unix and windows


=== "Solution"

    ```
    <option value="value_win" id="id_windows" > Windows </option>
    <option value="value_unix" id="id_unix" > UNIX </option>
    ```

You should now have a nearly finished dropdown box:

<div id="div_html" style="display:block" markdown="1">

```
<b>Please select the platform you wish to use for this exercise:<b>
<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
    <option value="value_win" id="id_windows" > Windows </option>
    <option value="value_unix" id="id_unix" > UNIX </option></select>
```

</div>

<div id="div_markdown" style="display:block" markdown="1">

```
**Please select the platform you wish to use for this exercise:**
<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
    <option value="value_win" id="id_windows" > Windows </option>
    <option value="value_unix" id="id_unix" > UNIX </option></select>
```

</div>

The last thing we need to do is set one of these options as the default selection.
This is different from the default page display we wrote into our function. The
function tells the website which parts of the page to display, what we need to
do here is set which option is shown as selected when the page loads. We can set
it to whatever we want, the page will always load with the default from the script,
but to make the page make sense to our users, the option it shows as selected,
should match the defaults we've set elsewhere. Since we've set Windows to the
defaults in the script, we should add 'selected' to the Windows option here:


<div id="div_html" style="display:block" markdown="1">


```
<b>Please select the platform you wish to use for this exercise:<b>
<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
    <option value="value_win" id="id_windows" selected> Windows </option>
    <option value="value_unix" id="id_unix" > UNIX </option></select>
```

</div>

<div id="div_markdown" style="display:block" markdown="1">

```
**Please select the platform you wish to use for this exercise:**
<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
    <option value="value_win" id="id_windows" selected> Windows </option>
    <option value="value_unix" id="id_unix" > UNIX </option></select>
```

</div>
