// See http://disqus.disqus.com/general_integration_help/
var disqus_developer = 1;
if (document.location.protocol == 'http:') {
    disqus_developer = 0;
}

/**
 * Give warnings about old documentation by detecting them in URL
 */
function summonDragons() {
    if(window.location.href.indexOf("/old/") >= 0) {
        $("#dragon-warning").show();
    }
}

/**
 * Make sure we get redirected to developer.ploen.org
 *
 */
function preferPrimaryDomain() {

    // Special cases
    if(window.location.protocol == "file:") {
        return;
    }

    if(window.location.hostname == "localhost") {
        return;
    }

    // Redirect rtd.org to developer.plone.org
    if(window.location.hostname != "developer.plone.org") {
        var href = window.location.href;
        href = href.replace(window.location.hostname, "developer.plone.org");
        href = href.replace("https://", "http://");
        window.location = href;
    }
}


preferPrimaryDomain();

// Some HTML/CSS transforms we can't do in pure CSS (waiting for CSS3 support)
$(document).ready(function () {

    // Style of .. verisonmodified:: , .. versionadded:: ...
    $('p:has(span.versionmodified)').addClass('p-versionmodified');

    // Style of the TOC tree
    $('ul:has(li.toctree-l1)').addClass('toctree');

    // Decorate links to other sites
    $('a[class~=external]').attr('target', '_blank');

    // Glossary terms
    $('em.xref').attr('title', _('See the definition'));

    // We need to have all <pre> in a <div class="highlight"> to mimic pygments
    // to work around a pygments (?) bug
    var suffixes = {cfg: null, python:null};
    for (var suffix in suffixes) {
        var pattern = 'div.highlight-' + suffix + ' > pre';
        $(pattern).wrap('<div class="highlight" />');
    }

    summonDragons();

    return true;
});

// Store editor pop-up help state in localStorage
// so it does not re-pop-up itself between page loads.
// Do not even to pretend to support IE gracefully.
(function($) {

    $(document).ready(function() {
        var box = $("#editor-trap");
        var klass = "toggled";
        var storageKey = "toggled";

        function toggle() {
            box.toggleClass(klass);
            // Store the toggle status in local storage as "has value string" or null
            window.localStorage.setItem(storageKey, box.hasClass(klass) ? "toggled" : "not-toggled");
        }

        box.click(toggle);

        // Check the persistent state of the editor pop-up
        // Note that localStorage does not necessarily support boolean values (ugh!)
        // http://stackoverflow.com/questions/3263161/cannot-set-boolean-values-in-localstorage
        var v = window.localStorage.getItem(storageKey);
        if(v == "toggled" || !v) {
          box.addClass(klass);
        }

    });

})(jQuery);
