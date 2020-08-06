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

window.onload = set_page_view_defaults;
</script>

**Please select the platform you wish to use for this exercise:**
<select id="id_platform" name="platformlist" onchange="change_content_by_platform('id_platform');return false;">
    <option value="value_mac" id="id_mac" > MacOS </option>
    <option value="value_win" id="id_windows" selected> Windows </option>
    <option value="value_unix" id="id_unix" > UNIX </option></select>


An introduction that is always displayed

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
