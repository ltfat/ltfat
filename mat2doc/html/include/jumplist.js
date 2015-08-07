$(document).ready(function(){
    // 
    // globalfun list and toolboxversion variables are defined
    if(typeof globalfunlist === 'undefined'){
        $('#jumplist').html('No function list found');
    }
    else{
        var reldots = "..";
        if(typeof rootdoc !== 'undefined'){
            reldots = ".";
        }


        // Key-sorting function
        var ksort = function ( src ) {
            var keys = Object.keys( src ),
                target = {};
            keys.sort();
            keys.forEach(function ( key ) {
                target[ key ] = src[ key ];
            });
            return target;
        };

        // Add endsWith to all strings 
        String.prototype.endsWith = function(suffix) {
            return this.indexOf(suffix, this.length - suffix.length) !== -1;
        };

        // Group by categories
        funlistByCat = {};
        for (var key in globalfunlist) {
            if (globalfunlist.hasOwnProperty(key)) {
                var tmp = globalfunlist[key];
                var pos = tmp.search('/');
                var cat = 'base'; 
                if(pos !== -1){
                    cat = tmp.substr(0,pos);
                }
                else{
                    tmp = cat + '/' + tmp;
                }

                cat = cat + ":";
                if(!funlistByCat.hasOwnProperty(cat)){
                    funlistByCat[cat] = {};
                }
                funlistByCat[cat][key]= tmp;
            }
        }

        funlistByCat = ksort(funlistByCat);

        //console.log(funlistByCat);
        // Build a combobox with function names
        var jstr = '';
        jstr += '<div class="input-group pull-right">';
        jstr += '<span class="input-group-addon"><b>Go to function</b></span>';
        //jstr += '<select class="form-control input-sm">';
        jstr += '<select class="selectpicker form-control input-sm" data-live-search="true">';
        jstr += '<option value="start/index">-----------------------</option>';

        // Get path from location, strip .html or _code.html
        var loc = window.location;
        var tocompare = ''; 
        var pos = loc.pathname.search('.html');
        var undcodesuf = '';
        if( pos !== -1){
            tocompare = loc.pathname.substr(0,pos);
            if(tocompare.endsWith('_code')){
                tocompare = tocompare.substr(0,tocompare.length-5);
                undcodesuf = '_code';
            }
        }


        for (var key in funlistByCat){
            // Add group label
            jstr += '<optgroup label="'+key+'">';
            tmpcatobj = ksort(funlistByCat[key]);
            // Add all elements, figure out which one is selected by comparing it to the
            // current path
            for(var fname in tmpcatobj){
                var selectedstr = '';
                var relpath = tmpcatobj[fname];
                if( tocompare.length>0 && tocompare.endsWith(relpath) ){
                    selectedstr = 'selected';
                }
                jstr += '<option value="'+relpath+'" '+selectedstr+'>'+fname+'</option>';
            }
            jstr += '</optgroup>';
        }

        jstr += '</select>';
        jstr += '</div>';
        $('#jumplist').html(jstr);


        // Redirect browser to new location on combobox change
        $('#jumplist select').change(function(){
            window.location.replace(reldots+'/' + this.value + undcodesuf + '.html');
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


