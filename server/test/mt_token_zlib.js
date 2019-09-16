var zlib = require("zlib");

var a = {"rId":100016,"ts":1555897671557,"cts":1555897996633,"brVD":[414,624],"brR":[[1242,1872],[1242,1872],24,24],"bI":["pages/index/index",""],"mT":[],"kT":[],"aT":[],"tT":[],"sign":"eJyVkEFuwyAURK9SecHOLmCD7QU36CIHiIQw/Dio2LiAnaSnL3GUqot2UbH684b5A4VWCUYfbjLdFhAYOX0WBGOCnFfmIRI0q82O2fg0LWoEaWcD1wz3IdpPEPROQpohiAYtTqWTD5MgNQrwIZOdchRjrOvbvueMcxRsfH+DDVxOic4aCDKCA52kUUmJY3EsUPQhPdeuqzWiJ5iA7nGpNYeyMXQoe8B1SSEf0xJgiqGLspOyudU4i1d0maTSaVVO5k42rQZEjWnb4I79ZH4eH5BQnIviFu90WTYI0fpZ1FVX7Tf0Xuhyzew+mu+3GtiS9y7u6m6yh7Of4aV9Obj1IT/TeMWr+q78WuqPNv/4gkluNtqU/ZR2pKGqLruB8bIZTrjsWDOUhA/M9Cfd1bwpvgBY6ayA"}
var b = JSON.stringify(a);

zlib.deflate(b, function(err, buffer) {
	if (!err) {
		console.log("deflate (%s): ", buffer.length, buffer.toString('base64'));
	}
});