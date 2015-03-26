$(document).ready(function(){
// 
// globalfun list and toolboxversion variables are defined
if(typeof globalfunlist === 'undefined'){
$('#jumplist').html('No function list found');
}
else{

// Build a combobox with function names
var jstr = '';
jstr += '<div class="input-group">';
jstr += '<span class="input-group-addon"><b>Go to function</b></span>';
jstr += '<form>';
jstr += '<select>';

for(var fname in globalfunlist){
    var relpath = globalfunlist[fname];
    jstr += '<option value="'+relpath+'">'+fname+'</option>';
}
jstr += '</select></form>';
jstr += '</div>';
$('#jumplist').html(jstr);


// Redirect browser to new location on combobox change
$('#jumplist select').change(function(){
    window.location.replace('../' + this.value + '.html');
});
}

// Add classes to the code-doc switch
$('#codeswitch a').addClass('btn btn-info btn-block');
// Add classes to seealso menu and move menutitle to the list
var menutitle = $('#seealso #menutitle').html();
$('#seealso #menutitle').html('<li><a href=""><b>'+menutitle +'</b></a></li>');
$('#seealso ul').addClass('nav doc-sidenav affix-top')
$('#seealso #menutitle > li').prependTo('#seealso ul');


});


