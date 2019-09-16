var now = require("performance-now");

function id() {
	var e = (new Date).getTime();
	return (e += now()), "xxxxxxxxxxxx4xxxyxxxxxxxxxxxxxxx".replace(/[xy]/g, function(t) {
		var r = (e + 16 * Math.random()) % 16 | 0;
		return e = Math.floor(e / 16), ("x" === t ? r : 3 & r | 8).toString(16)
	}).toUpperCase() + "|" + Date.now()
};

function getId() {
	var t = "LoginService.loginByUsername"
	var r = {
		captchaCode: "bxpi",
		loginedSessionIds: [],
		password: '11111',
		username: '22222'
	}
	m = function(e, t) {
		return true
	}
	S = function() {
		return function(e, t) {
			if (Array.isArray(e))
				return e;
			if (Symbol.iterator in Object(e))
				return function(e, t) {
					var r = [],
						n = !0,
						i = !1,
						o = void 0;
					try {
						for (var a, s = e[Symbol.iterator](); !(n = (a = s.next()).done) && (r.push(a.value), !t || r.length !== t); n = !0)
							;
					}
					catch (e) {
						i = !0,
							o = e
					}
					finally {
						try {
							!n && s.return && s.return()
						}
						finally {
							if (i)
								throw o
						}
					}
					return r
				}(e, t);
			throw new TypeError("Invalid attempt to destructure non-iterable instance")
		}
	}()

	function e(t, r) {
		var n = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : {};
		m(this, e);
		var i = t.split("."),
			o = S(i, 2),
			a = o[0],
			s = o[1]
		this.id = function() {
			var e = new Date().getTime();
			return true && "function" == "function" && (e += 2637349.000000002),
			"xxxxxxxxxxxx4xxxyxxxxxxxxxxxxxxx".replace(/[xy]/g, function(t) {
				var r = (e + 16 * Math.random()) % 16 | 0;
				return e = Math.floor(e / 16),
					("x" === t ? r : 3 & r | 8).toString(16)
			}).toUpperCase() + "|" + Date.now()
		}(),
			this.metas = n.metas,
			this.service = a,
			this.method = s,
			this.params = r,
			this.ncp = n.ncp
		return this.id
	}
	var self = e(t, r)
	return self
}

console.log(id())
console.log(getId())