$(document).ready(function() {
	$(".sidebar").on("activate", function(){
		$(".usermenu li a i").removeClass("icon-white");
		$(".usermenu li.active a i").addClass("icon-white")
	});

	//fancy scrolling animation	
	$('.usermenu li a').on('click', function(e) {
	   // prevent default anchor click behavior
	   e.preventDefault();

	   // store hash
	   var hash = this.hash;

	   // animate
	   $('html, body').animate({
	       scrollTop: $(this.hash).offset().top
	     }, 500, function(){
	       // when done, add hash to url
	       // (default click behaviour)
	       window.location.hash = hash;
	     });
	});	

	//see https://github.com/twitter/bootstrap/issues/6350
	$('[data-clampedwidth]').each(function () {
		var elem = $(this);
		var parentPanel = elem.data('clampedwidth');
		var resizeFn = function () {
			var sideBarNavWidth = $(parentPanel).width() - parseInt(elem.css('paddingLeft')) - parseInt(elem.css('paddingRight')) - parseInt(elem.css('marginLeft')) - parseInt(elem.css('marginRight')) - parseInt(elem.css('borderLeftWidth')) - parseInt(elem.css('borderRightWidth'));
			elem.css('width', sideBarNavWidth);
		};

		resizeFn();
		$(window).resize(resizeFn);
	});	

	$("#modallink").click(function(e){
		if(window.innerWidth > 768) {
			e.preventDefault();
			$('#twittermodal').modal('toggle');
			return false;
		}
	});
});
