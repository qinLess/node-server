var babelHelpers = {}

var ts = Date.now();
var cts = Date.now() + 13254;

var window = {
    seed: {
        env: "pro"
    },
    innerWidth: 1920,
    innerHeight: 1080,
    location: {
        href: ""
    },
    bi: null,
    Yoda: {},
    behavior: '',
    token: '',
    ts: ts,
    cts: cts,
    token_btoa: '',
    img_token: '',
    Navigator: {},
    Image: {},
    Screen: {},
    Audio: {},
    Location: {},
    sleep: function (numberMillis) {
        var now = new Date();
        var exitTime = now.getTime() + numberMillis;
        while (true) {
            now = new Date();
            if (now.getTime() > exitTime)
                return;
        }
    },
    btoa: function (str) {
        var buffer;
        if (str instanceof Buffer) {
            buffer = str;
        } else {
            buffer = Buffer.from(str.toString(), 'binary');
        }
        return buffer.toString('base64');
    },
    atob: function (str) {
        return Buffer.from(str, 'base64').toString('binary');
    },

    point: function() {
        // 滑块
        var a = 16 * 2.735
        // 滑框
        var b = 414 * 0.9
        // 边距
        var c = (414 - b) / 2
        // 最大值
        var d = b + c
        // 最小值1
        var e = b - a + c
        // y 轴 最大值
        var f = 736

        var trace = [];

        var F = 0

        while (true) {
            // y 轴 上下浮动
            var g = Math.ceil(Math.random() * a);
            if (F > d) {
                break
            }
            F = F + ((a / 2 + c) - Math.ceil(Math.random() * 10))
            var S = f + g;
            var T = Math.floor(Math.random() * (20000 - 15000) + 1) + 15000;
            trace.push([0, parseInt(F), parseInt(S) / 2, T]);
        }

        return trace
    },

    mT: function () {
        var mt = [];
        // var k = Math.floor(Math.random() * ((60 - 40) + 1)) + 40;
        var k = Math.floor(Math.random() * 60);
        var arr = ['BUTTON', 'INPUT', 'DIV'];
        for (var i = 1; i <= k; i++) {
            window.sleep(0.5)
            var page_x = Math.floor(Math.random() * ((1283 - 482) + 1)) + 482;
            var page_y = Math.floor(Math.random() * ((295 - 245) + 1)) + 245;
            var t = Date.now() - window.ts;
            // 			var s = arr[Math.floor((Math.random() * arr.length))];
            mt.push(`${page_x}, ${page_y}, ${t}`)
        }
        mt.length > 60 && (mt = mt.slice(0, 60))
        return mt
    },
    kT: function () {
        var k_T = []
        for (var t = 1; t <= 5; t++) {
            window.sleep(0.2);
            var i = Date.now() - window.ts;
            var s = Math.floor(Math.random() * ((90) - (65 + 1))) + 65;
            var r = String(s).charCodeAt();
            k_T.push(`${r}, INPUT, ${i}`)

        }
        return k_T
    },
    aT: function () {
        var a_T = [];
        var list = ['BUTTON', 'INPUT', 'DIV', 'A', 'HTML'];
        for (var t = 1; t < 10; t++) {
            var page_x = Math.floor(Math.random() * ((1283) - (482 + 1))) + 482;
            var page_y = Math.floor(Math.random() * ((295) - (245 + 1))) + 245;
            var i = Date.now() - window.ts;
            var s = list[Math.floor(Math.random() * list.length)];
            a_T.push(`${page_x}, ${page_y}, ${s}, ${i}`)
        }
        return a_T
    },
    tT: function () {
        var t_T = [];
        var k = Math.floor(Math.random() * ((20 - 15) + 1)) + 15;
        for (var t = 1; t < k; t++) {
            var page_x = Math.floor(Math.random() * ((360) - (60 + 1))) + 60;
            var page_y = Math.floor(Math.random() * ((370) - (350 + 1))) + 350;
            var i = Math.floor(Math.random() * ((18000) - (17000 + 1))) + 17000;
            t_T.push(`${page_x}, ${page_y}, 1, ${i}`)
        }
        return t_T
    }
};

var navigator = {
    userAgent: ''
};

var document = {
    documentElement: {
        clientWidth: 800,
        clientHeight: 600
    }
};
var screen = {
    width: 1920,
    height: 1080,
    availWidth: 1858,
    availHeight: 1080,
    colorDepth: 24,
    pixelDepth: 24
};


function abcdefg(user_agent) {
    navigator.userAgent = user_agent
    window.navigator = navigator
    window.document = document
    window.screen = screen

    ;(function() {
            var g;
            if (typeof window !== "undefined") {
                g = window
            } else if (typeof global !== "undefined") {
                g = global
            } else if (typeof self !== "undefined") {
                g = self
            } else {
                g = this
            }
            g.babelHelpers = babelHelpers;
            babelHelpers.typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function(obj) {
                    return typeof obj
                }
                : function(obj) {
                    return obj && typeof Symbol === "function" && obj.constructor === Symbol ? "symbol" : typeof obj
                }
            ;
            babelHelpers.jsx = function() {
                var REACT_ELEMENT_TYPE = typeof Symbol === "function" && Symbol.for && Symbol.for("react.element") || 60103;
                return function createRawReactElement(type, props, key, children) {
                    var defaultProps = type && type.defaultProps;
                    var childrenLength = arguments.length - 3;
                    if (!props && childrenLength !== 0) {
                        props = {}
                    }
                    if (props && defaultProps) {
                        for (var propName in defaultProps) {
                            if (props[propName] === void 0) {
                                props[propName] = defaultProps[propName]
                            }
                        }
                    } else if (!props) {
                        props = defaultProps || {}
                    }
                    if (childrenLength === 1) {
                        props.children = children
                    } else if (childrenLength > 1) {
                        var childArray = Array(childrenLength);
                        for (var i = 0; i < childrenLength; i++) {
                            childArray[i] = arguments[i + 3]
                        }
                        props.children = childArray
                    }
                    return {
                        $$typeof: REACT_ELEMENT_TYPE,
                        type: type,
                        key: key === undefined ? null : "" + key,
                        ref: null,
                        props: props,
                        _owner: null
                    }
                }
            }();
            babelHelpers.asyncToGenerator = function(fn) {
                return function() {
                    var gen = fn.apply(this, arguments);
                    return new Promise(function(resolve, reject) {
                            function step(key, arg) {
                                try {
                                    var info = gen[key](arg);
                                    var value = info.value
                                } catch (error) {
                                    reject(error);
                                    return
                                }
                                if (info.done) {
                                    resolve(value)
                                } else {
                                    return Promise.resolve(value).then(function(value) {
                                        return step("next", value)
                                    }, function(err) {
                                        return step("throw", err)
                                    })
                                }
                            }
                            return step("next")
                        }
                    )
                }
            }
            ;
            babelHelpers.classCallCheck = function(instance, Constructor) {
                if (!(instance instanceof Constructor)) {
                    throw new TypeError("Cannot call a class as a function")
                }
            }
            ;
            babelHelpers.createClass = function() {
                function defineProperties(target, props) {
                    for (var i = 0; i < props.length; i++) {
                        var descriptor = props[i];
                        descriptor.enumerable = descriptor.enumerable || false;
                        descriptor.configurable = true;
                        if ("value"in descriptor)
                            descriptor.writable = true;
                        Object.defineProperty(target, descriptor.key, descriptor)
                    }
                }
                return function(Constructor, protoProps, staticProps) {
                    if (protoProps)
                        defineProperties(Constructor.prototype, protoProps);
                    if (staticProps)
                        defineProperties(Constructor, staticProps);
                    return Constructor
                }
            }();
            babelHelpers.defineEnumerableProperties = function(obj, descs) {
                for (var key in descs) {
                    var desc = descs[key];
                    desc.configurable = desc.enumerable = true;
                    if ("value"in desc)
                        desc.writable = true;
                    Object.defineProperty(obj, key, desc)
                }
                return obj
            }
            ;
            babelHelpers.defaults = function(obj, defaults) {
                var keys = Object.getOwnPropertyNames(defaults);
                for (var i = 0; i < keys.length; i++) {
                    var key = keys[i];
                    var value = Object.getOwnPropertyDescriptor(defaults, key);
                    if (value && value.configurable && obj[key] === undefined) {
                        Object.defineProperty(obj, key, value)
                    }
                }
                return obj
            }
            ;
            babelHelpers.defineProperty = function(obj, key, value) {
                if (key in obj) {
                    Object.defineProperty(obj, key, {
                        value: value,
                        enumerable: true,
                        configurable: true,
                        writable: true
                    })
                } else {
                    obj[key] = value
                }
                return obj
            }
            ;
            babelHelpers.extends = Object.assign || function(target) {
                for (var i = 1; i < arguments.length; i++) {
                    var source = arguments[i];
                    for (var key in source) {
                        if (Object.prototype.hasOwnProperty.call(source, key)) {
                            target[key] = source[key]
                        }
                    }
                }
                return target
            }
            ;
            babelHelpers.get = function get(object, property, receiver) {
                if (object === null)
                    object = Function.prototype;
                var desc = Object.getOwnPropertyDescriptor(object, property);
                if (desc === undefined) {
                    var parent = Object.getPrototypeOf(object);
                    if (parent === null) {
                        return undefined
                    } else {
                        return get(parent, property, receiver)
                    }
                } else if ("value"in desc) {
                    return desc.value
                } else {
                    var getter = desc.get;
                    if (getter === undefined) {
                        return undefined
                    }
                    return getter.call(receiver)
                }
            }
            ;
            babelHelpers.inherits = function(subClass, superClass) {
                if (typeof superClass !== "function" && superClass !== null) {
                    throw new TypeError("Super expression must either be null or a function, not " + typeof superClass)
                }
                subClass.prototype = Object.create(superClass && superClass.prototype, {
                    constructor: {
                        value: subClass,
                        enumerable: false,
                        writable: true,
                        configurable: true
                    }
                });
                if (superClass)
                    Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass
            }
            ;
            babelHelpers.instanceof = function(left, right) {
                if (right != null && typeof Symbol !== "undefined" && right[Symbol.hasInstance]) {
                    return right[Symbol.hasInstance](left)
                } else {
                    return left instanceof right
                }
            }
            ;
            babelHelpers.interopRequireDefault = function(obj) {
                return obj && obj.__esModule ? obj : {
                    "default": obj
                }
            }
            ;
            babelHelpers.interopRequireWildcard = function(obj) {
                if (obj && obj.__esModule) {
                    return obj
                } else {
                    var newObj = {};
                    if (obj != null) {
                        for (var key in obj) {
                            if (Object.prototype.hasOwnProperty.call(obj, key))
                                newObj[key] = obj[key]
                        }
                    }
                    newObj.default = obj;
                    return newObj
                }
            }
            ;
            babelHelpers.newArrowCheck = function(innerThis, boundThis) {
                if (innerThis !== boundThis) {
                    throw new TypeError("Cannot instantiate an arrow function")
                }
            }
            ;
            babelHelpers.objectDestructuringEmpty = function(obj) {
                if (obj == null)
                    throw new TypeError("Cannot destructure undefined")
            }
            ;
            babelHelpers.objectWithoutProperties = function(obj, keys) {
                var target = {};
                for (var i in obj) {
                    if (keys.indexOf(i) >= 0)
                        continue;
                    if (!Object.prototype.hasOwnProperty.call(obj, i))
                        continue;
                    target[i] = obj[i]
                }
                return target
            }
            ;
            babelHelpers.possibleConstructorReturn = function(self, call) {
                if (!self) {
                    throw new ReferenceError("this hasn't been initialised - super() hasn't been called")
                }
                return call && (typeof call === "object" || typeof call === "function") ? call : self
            }
            ;
            babelHelpers.selfGlobal = typeof global === "undefined" ? self : global;
            babelHelpers.set = function set(object, property, value, receiver) {
                var desc = Object.getOwnPropertyDescriptor(object, property);
                if (desc === undefined) {
                    var parent = Object.getPrototypeOf(object);
                    if (parent !== null) {
                        set(parent, property, value, receiver)
                    }
                } else if ("value"in desc && desc.writable) {
                    desc.value = value
                } else {
                    var setter = desc.set;
                    if (setter !== undefined) {
                        setter.call(receiver, value)
                    }
                }
                return value
            }
            ;
            babelHelpers.slicedToArray = function() {
                function sliceIterator(arr, i) {
                    var _arr = [];
                    var _n = true;
                    var _d = false;
                    var _e = undefined;
                    try {
                        for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) {
                            _arr.push(_s.value);
                            if (i && _arr.length === i)
                                break
                        }
                    } catch (err) {
                        _d = true;
                        _e = err
                    } finally {
                        try {
                            if (!_n && _i["return"])
                                _i["return"]()
                        } finally {
                            if (_d)
                                throw _e
                        }
                    }
                    return _arr
                }
                return function(arr, i) {
                    if (Array.isArray(arr)) {
                        return arr
                    } else if (Symbol.iterator in Object(arr)) {
                        return sliceIterator(arr, i)
                    } else {
                        throw new TypeError("Invalid attempt to destructure non-iterable instance")
                    }
                }
            }();
            babelHelpers.slicedToArrayLoose = function(arr, i) {
                if (Array.isArray(arr)) {
                    return arr
                } else if (Symbol.iterator in Object(arr)) {
                    var _arr = [];
                    for (var _iterator = arr[Symbol.iterator](), _step; !(_step = _iterator.next()).done; ) {
                        _arr.push(_step.value);
                        if (i && _arr.length === i)
                            break
                    }
                    return _arr
                } else {
                    throw new TypeError("Invalid attempt to destructure non-iterable instance")
                }
            }
            ;
            babelHelpers.taggedTemplateLiteral = function(strings, raw) {
                return Object.freeze(Object.defineProperties(strings, {
                    raw: {
                        value: Object.freeze(raw)
                    }
                }))
            }
            ;
            babelHelpers.taggedTemplateLiteralLoose = function(strings, raw) {
                strings.raw = raw;
                return strings
            }
            ;
            babelHelpers.temporalRef = function(val, name, undef) {
                if (val === undef) {
                    throw new ReferenceError(name + " is not defined - temporal dead zone")
                } else {
                    return val
                }
            }
            ;
            babelHelpers.temporalUndefined = {};
            babelHelpers.toArray = function(arr) {
                return Array.isArray(arr) ? arr : Array.from(arr)
            }
            ;
            babelHelpers.toConsumableArray = function(arr) {
                if (Array.isArray(arr)) {
                    for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++)
                        arr2[i] = arr[i];
                    return arr2
                } else {
                    return Array.from(arr)
                }
            }
        }
    )();
    /* Yoda loader for mobile | 2019-3-7 19:02:55 */
    !function(t) {
        function e(r) {
            if (n[r])
                return n[r].exports;
            var i = n[r] = {
                exports: {},
                id: r,
                loaded: !1
            };
            return t[r].call(i.exports, i, i.exports, e),
                i.loaded = !0,
                i.exports
        }
        var n = {};
        return e.m = t,
            e.c = n,
            e.p = "",
            e(0)
    }([function(t, e, n) {
        "use strict";
        n(1)
    }, function(t, e, n) {
            "use strict";
            n(2);
            var r = n(3)
                , i = babelHelpers.interopRequireDefault(r)
                , o = n(4)
                , a = babelHelpers.interopRequireDefault(o)
                , s = n(5)
                , u = babelHelpers.interopRequireDefault(s)
                , l = n(6)
                , c = babelHelpers.interopRequireDefault(l)
                , f = n(30)
                , h = babelHelpers.interopRequireDefault(f)
                , d = n(35)
                , p = babelHelpers.interopRequireDefault(d);
            (0,
                a.default)(window),
                window.Yoda.request = u.default,
                window.Yoda.Promise = c.default,
                window.Yoda.Ballade = h.default,
                window.Yoda.Adapter = p.default,
                window.Yoda.encode = i.default.encode,
                window.Yoda.decode = i.default.decode
        }
        , function(t, e) {
            "use strict";
            window.YODA_CONFIG = {},
                window.YODA_CONFIG.__APP_NAME__ = "yoda",
                window.YODA_CONFIG.__API_URL__ = "https://verify.meituan.com",
                window.YODA_CONFIG.__ENV__ = "production"
        }
        , function(t, e) {
            "use strict";
            var n = {};
            n.PADCHAR = "=",
                n.ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/",
                n.makeDOMException = function() {
                    try {
                        return new DOMException(DOMException.INVALID_CHARACTER_ERR)
                    } catch (e) {
                        var t = new Error("DOM Exception 5");
                        return t.code = t.number = 5,
                            t.name = t.description = "INVALID_CHARACTER_ERR",
                            t.toString = function() {
                                return "Error: " + t.name + ": " + t.message
                            }
                            ,
                            t
                    }
                }
                ,
                n.getbyte64 = function(t, e) {
                    var r = n.ALPHA.indexOf(t.charAt(e));
                    if (r === -1)
                        throw n.makeDOMException();
                    return r
                }
                ,
                n.decode = function(t) {
                    t = "" + t;
                    var e, r, i, o = n.getbyte64, a = t.length;
                    if (0 === a)
                        return t;
                    if (a % 4 !== 0)
                        throw n.makeDOMException();
                    e = 0,
                    t.charAt(a - 1) === n.PADCHAR && (e = 1,
                    t.charAt(a - 2) === n.PADCHAR && (e = 2),
                        a -= 4);
                    var s = [];
                    for (r = 0; r < a; r += 4)
                        i = o(t, r) << 18 | o(t, r + 1) << 12 | o(t, r + 2) << 6 | o(t, r + 3),
                            s.push(String.fromCharCode(i >> 16, i >> 8 & 255, 255 & i));
                    switch (e) {
                        case 1:
                            i = o(t, r) << 18 | o(t, r + 1) << 12 | o(t, r + 2) << 6,
                                s.push(String.fromCharCode(i >> 16, i >> 8 & 255));
                            break;
                        case 2:
                            i = o(t, r) << 18 | o(t, r + 1) << 12,
                                s.push(String.fromCharCode(i >> 16))
                    }
                    return s.join("")
                }
                ,
                n.getbyte = function(t, e) {
                    var r = t.charCodeAt(e);
                    if (r > 255)
                        throw n.makeDOMException();
                    return r
                }
                ,
                n.encode = function(t) {
                    if (1 !== arguments.length)
                        throw new SyntaxError("Not enough arguments");
                    var e, r, i = n.PADCHAR, o = n.ALPHA, a = n.getbyte, s = [];
                    t = "" + t;
                    var u = t.length - t.length % 3;
                    if (0 === t.length)
                        return t;
                    for (e = 0; e < u; e += 3)
                        r = a(t, e) << 16 | a(t, e + 1) << 8 | a(t, e + 2),
                            s.push(o.charAt(r >> 18)),
                            s.push(o.charAt(r >> 12 & 63)),
                            s.push(o.charAt(r >> 6 & 63)),
                            s.push(o.charAt(63 & r));
                    switch (t.length - u) {
                        case 1:
                            r = a(t, e) << 16,
                                s.push(o.charAt(r >> 18) + o.charAt(r >> 12 & 63) + i + i);
                            break;
                        case 2:
                            r = a(t, e) << 16 | a(t, e + 1) << 8,
                                s.push(o.charAt(r >> 18) + o.charAt(r >> 12 & 63) + o.charAt(r >> 6 & 63) + i)
                    }
                    return s.join("")
                }
                ,
            window.btoa || (window.btoa = n.encode),
            window.atob || (window.atob = n.decode),
                t.exports = n
        }
        , function(t, e) {
            "use strict";
            Object.defineProperty(e, "__esModule", {
                value: !0
            }),
                e.default = function(t) {
                    function e(t) {
                        switch ("undefined" == typeof t ? "undefined" : babelHelpers.typeof(t)) {
                            case "undefined":
                                return "undefined";
                            case "boolean":
                                return "boolean";
                            case "number":
                                return "number";
                            case "string":
                                return "string";
                            default:
                                return null === t ? "null" : "object"
                        }
                    }
                    function n(t) {
                        return Object.prototype.toString.call(t).replace(/^\[object *|\]$/g, "")
                    }
                    function r(t) {
                        return "function" == typeof t
                    }
                    function i(t) {
                        if (null === t || t === P)
                            throw TypeError();
                        return Object(t)
                    }
                    function o(t) {
                        return t >> 0
                    }
                    function a(t) {
                        return t >>> 0
                    }
                    function s(e) {
                        function n(t) {
                            Object.defineProperty(e, t, {
                                get: function() {
                                    return e._getter(t)
                                },
                                set: function(n) {
                                    e._setter(t, n)
                                },
                                enumerable: !0,
                                configurable: !1
                            })
                        }
                        if (!("TYPED_ARRAY_POLYFILL_NO_ARRAY_ACCESSORS"in t)) {
                            if (e.length > C)
                                throw RangeError("Array too large for polyfill");
                            var r;
                            for (r = 0; r < e.length; r += 1)
                                n(r)
                        }
                    }
                    function u(t, e) {
                        var n = 32 - e;
                        return t << n >> n
                    }
                    function l(t, e) {
                        var n = 32 - e;
                        return t << n >>> n
                    }
                    function c(t) {
                        return [255 & t]
                    }
                    function f(t) {
                        return u(t[0], 8)
                    }
                    function h(t) {
                        return [255 & t]
                    }
                    function d(t) {
                        return l(t[0], 8)
                    }
                    function p(t) {
                        return t = M(Number(t)),
                            [t < 0 ? 0 : t > 255 ? 255 : 255 & t]
                    }
                    function _(t) {
                        return [255 & t, t >> 8 & 255]
                    }
                    function y(t) {
                        return u(t[1] << 8 | t[0], 16)
                    }
                    function b(t) {
                        return [255 & t, t >> 8 & 255]
                    }
                    function g(t) {
                        return l(t[1] << 8 | t[0], 16)
                    }
                    function v(t) {
                        return [255 & t, t >> 8 & 255, t >> 16 & 255, t >> 24 & 255]
                    }
                    function w(t) {
                        return u(t[3] << 24 | t[2] << 16 | t[1] << 8 | t[0], 32)
                    }
                    function m(t) {
                        return [255 & t, t >> 8 & 255, t >> 16 & 255, t >> 24 & 255]
                    }
                    function E(t) {
                        return l(t[3] << 24 | t[2] << 16 | t[1] << 8 | t[0], 32)
                    }
                    function k(t, e, n) {
                        function r(t) {
                            var e = N(t)
                                , n = t - e;
                            return n < .5 ? e : n > .5 ? e + 1 : e % 2 ? e + 1 : e
                        }
                        var i, o, a, s = (1 << e - 1) - 1;
                        if (t !== t)
                            o = (1 << e) - 1,
                                a = B(2, n - 1),
                                i = 0;
                        else if (t === 1 / 0 || t === -(1 / 0))
                            o = (1 << e) - 1,
                                a = 0,
                                i = t < 0 ? 1 : 0;
                        else if (0 === t)
                            o = 0,
                                a = 0,
                                i = 1 / t === -(1 / 0) ? 1 : 0;
                        else if (i = t < 0,
                            t = z(t),
                        t >= B(2, 1 - s)) {
                            o = L(N(R(t) / S), 1023);
                            var u = t / B(2, o);
                            u < 1 && (o -= 1,
                                u *= 2),
                            u >= 2 && (o += 1,
                                u /= 2);
                            var l = B(2, n);
                            a = r(u * l) - l,
                                o += s,
                            a / l >= 1 && (o += 1,
                                a = 0),
                            o > 2 * s && (o = (1 << e) - 1,
                                a = 0)
                        } else
                            o = 0,
                                a = r(t / B(2, 1 - s - n));
                        var c, f = [];
                        for (c = n; c; c -= 1)
                            f.push(a % 2 ? 1 : 0),
                                a = N(a / 2);
                        for (c = e; c; c -= 1)
                            f.push(o % 2 ? 1 : 0),
                                o = N(o / 2);
                        f.push(i ? 1 : 0),
                            f.reverse();
                        for (var h = f.join(""), d = []; h.length; )
                            d.unshift(parseInt(h.substring(0, 8), 2)),
                                h = h.substring(8);
                        return d
                    }
                    function O(t, e, n) {
                        var r, i, o, a, s, u, l, c, f = [];
                        for (r = 0; r < t.length; ++r)
                            for (o = t[r],
                                    i = 8; i; i -= 1)
                                f.push(o % 2 ? 1 : 0),
                                    o >>= 1;
                        return f.reverse(),
                            a = f.join(""),
                            s = (1 << e - 1) - 1,
                            u = parseInt(a.substring(0, 1), 2) ? -1 : 1,
                            l = parseInt(a.substring(1, 1 + e), 2),
                            c = parseInt(a.substring(1 + e), 2),
                            l === (1 << e) - 1 ? 0 !== c ? NaN : u * (1 / 0) : l > 0 ? u * B(2, l - s) * (1 + c / B(2, n)) : 0 !== c ? u * B(2, -(s - 1)) * (c / B(2, n)) : u < 0 ? -0 : 0
                    }
                    function T(t) {
                        return O(t, 11, 52)
                    }
                    function j(t) {
                        return k(t, 11, 52)
                    }
                    function A(t) {
                        return O(t, 8, 23)
                    }
                    function x(t) {
                        return k(t, 8, 23)
                    }
                    var P = void 0
                        , C = 1e5
                        , S = Math.LN2
                        , z = Math.abs
                        , N = Math.floor
                        , R = Math.log
                        , I = Math.max
                        , L = Math.min
                        , B = Math.pow
                        , M = Math.round;
                    !function() {
                        var t = Object.defineProperty
                            , e = !function() {
                            try {
                                return Object.defineProperty({}, "x", {})
                            } catch (t) {
                                return !1
                            }
                        }();
                        t && !e || (Object.defineProperty = function(e, n, r) {
                                if (t)
                                    try {
                                        return t(e, n, r)
                                    } catch (t) {}
                                if (e !== Object(e))
                                    throw TypeError("Object.defineProperty called on non-object");
                                return Object.prototype.__defineGetter__ && "get"in r && Object.prototype.__defineGetter__.call(e, n, r.get),
                                Object.prototype.__defineSetter__ && "set"in r && Object.prototype.__defineSetter__.call(e, n, r.set),
                                "value"in r && (e[n] = r.value),
                                    e
                            }
                        )
                    }(),
                        function() {
                            function u(t) {
                                if (t = o(t),
                                t < 0)
                                    throw RangeError("ArrayBuffer size is not a small enough positive integer.");
                                Object.defineProperty(this, "byteLength", {
                                    value: t
                                }),
                                    Object.defineProperty(this, "_bytes", {
                                        value: Array(t)
                                    });
                                for (var e = 0; e < t; e += 1)
                                    this._bytes[e] = 0
                            }
                            function l() {
                                if (!arguments.length || "object" !== babelHelpers.typeof(arguments[0]))
                                    return function(t) {
                                        if (t = o(t),
                                        t < 0)
                                            throw RangeError("length is not a small enough positive integer.");
                                        Object.defineProperty(this, "length", {
                                            value: t
                                        }),
                                            Object.defineProperty(this, "byteLength", {
                                                value: t * this.BYTES_PER_ELEMENT
                                            }),
                                            Object.defineProperty(this, "buffer", {
                                                value: new u(this.byteLength)
                                            }),
                                            Object.defineProperty(this, "byteOffset", {
                                                value: 0
                                            })
                                    }
                                        .apply(this, arguments);
                                if (arguments.length >= 1 && "object" === e(arguments[0]) && arguments[0]instanceof l)
                                    return function(t) {
                                        if (this.constructor !== t.constructor)
                                            throw TypeError();
                                        var e = t.length * this.BYTES_PER_ELEMENT;
                                        Object.defineProperty(this, "buffer", {
                                            value: new u(e)
                                        }),
                                            Object.defineProperty(this, "byteLength", {
                                                value: e
                                            }),
                                            Object.defineProperty(this, "byteOffset", {
                                                value: 0
                                            }),
                                            Object.defineProperty(this, "length", {
                                                value: t.length
                                            });
                                        for (var n = 0; n < this.length; n += 1)
                                            this._setter(n, t._getter(n))
                                    }
                                        .apply(this, arguments);
                                if (arguments.length >= 1 && "object" === e(arguments[0]) && !(arguments[0]instanceof l) && !(arguments[0]instanceof u || "ArrayBuffer" === n(arguments[0])))
                                    return function(t) {
                                        var e = t.length * this.BYTES_PER_ELEMENT;
                                        Object.defineProperty(this, "buffer", {
                                            value: new u(e)
                                        }),
                                            Object.defineProperty(this, "byteLength", {
                                                value: e
                                            }),
                                            Object.defineProperty(this, "byteOffset", {
                                                value: 0
                                            }),
                                            Object.defineProperty(this, "length", {
                                                value: t.length
                                            });
                                        for (var n = 0; n < this.length; n += 1) {
                                            var r = t[n];
                                            this._setter(n, Number(r))
                                        }
                                    }
                                        .apply(this, arguments);
                                if (arguments.length >= 1 && "object" === e(arguments[0]) && (arguments[0]instanceof u || "ArrayBuffer" === n(arguments[0])))
                                    return function(t, e, n) {
                                        if (e = a(e),
                                        e > t.byteLength)
                                            throw RangeError("byteOffset out of range");
                                        if (e % this.BYTES_PER_ELEMENT)
                                            throw RangeError("buffer length minus the byteOffset is not a multiple of the element size.");
                                        if (n === P) {
                                            var r = t.byteLength - e;
                                            if (r % this.BYTES_PER_ELEMENT)
                                                throw RangeError("length of buffer minus byteOffset not a multiple of the element size");
                                            n = r / this.BYTES_PER_ELEMENT
                                        } else
                                            n = a(n),
                                                r = n * this.BYTES_PER_ELEMENT;
                                        if (e + r > t.byteLength)
                                            throw RangeError("byteOffset and length reference an area beyond the end of the buffer");
                                        Object.defineProperty(this, "buffer", {
                                            value: t
                                        }),
                                            Object.defineProperty(this, "byteLength", {
                                                value: r
                                            }),
                                            Object.defineProperty(this, "byteOffset", {
                                                value: e
                                            }),
                                            Object.defineProperty(this, "length", {
                                                value: n
                                            })
                                    }
                                        .apply(this, arguments);
                                throw TypeError()
                            }
                            function k(t, e, n) {
                                var r = function t() {
                                    Object.defineProperty(this, "constructor", {
                                        value: t
                                    }),
                                        l.apply(this, arguments),
                                        s(this)
                                };
                                "__proto__"in r ? r.__proto__ = l : (r.from = l.from,
                                    r.of = l.of),
                                    r.BYTES_PER_ELEMENT = t;
                                var i = function() {};
                                return i.prototype = O,
                                    r.prototype = new i,
                                    Object.defineProperty(r.prototype, "BYTES_PER_ELEMENT", {
                                        value: t
                                    }),
                                    Object.defineProperty(r.prototype, "_pack", {
                                        value: e
                                    }),
                                    Object.defineProperty(r.prototype, "_unpack", {
                                        value: n
                                    }),
                                    r
                            }
                            t.ArrayBuffer = t.ArrayBuffer || u,
                                Object.defineProperty(l, "from", {
                                    value: function(t) {
                                        return new this(t)
                                    }
                                }),
                                Object.defineProperty(l, "of", {
                                    value: function() {
                                        return new this(arguments)
                                    }
                                });
                            var O = {};
                            l.prototype = O,
                                Object.defineProperty(l.prototype, "_getter", {
                                    value: function(t) {
                                        if (arguments.length < 1)
                                            throw SyntaxError("Not enough arguments");
                                        if (t = a(t),
                                        t >= this.length)
                                            return P;
                                        var e, n, r = [];
                                        for (e = 0,
                                                n = this.byteOffset + t * this.BYTES_PER_ELEMENT; e < this.BYTES_PER_ELEMENT; e += 1,
                                                n += 1)
                                            r.push(this.buffer._bytes[n]);
                                        return this._unpack(r)
                                    }
                                }),
                                Object.defineProperty(l.prototype, "get", {
                                    value: l.prototype._getter
                                }),
                                Object.defineProperty(l.prototype, "_setter", {
                                    value: function(t, e) {
                                        if (arguments.length < 2)
                                            throw SyntaxError("Not enough arguments");
                                        if (t = a(t),
                                            !(t >= this.length)) {
                                            var n, r, i = this._pack(e);
                                            for (n = 0,
                                                    r = this.byteOffset + t * this.BYTES_PER_ELEMENT; n < this.BYTES_PER_ELEMENT; n += 1,
                                                    r += 1)
                                                this.buffer._bytes[r] = i[n]
                                        }
                                    }
                                }),
                                Object.defineProperty(l.prototype, "constructor", {
                                    value: l
                                }),
                                Object.defineProperty(l.prototype, "copyWithin", {
                                    value: function(t, e) {
                                        var n = arguments[2]
                                            , r = i(this)
                                            , s = r.length
                                            , u = a(s);
                                        u = I(u, 0);
                                        var l, c = o(t);
                                        l = c < 0 ? I(u + c, 0) : L(c, u);
                                        var f, h = o(e);
                                        f = h < 0 ? I(u + h, 0) : L(h, u);
                                        var d;
                                        d = n === P ? u : o(n);
                                        var p;
                                        p = d < 0 ? I(u + d, 0) : L(d, u);
                                        var _, y = L(p - f, u - l);
                                        for (f < l && l < f + y ? (_ = -1,
                                            f = f + y - 1,
                                            l = l + y - 1) : _ = 1; y > 0; )
                                            r._setter(l, r._getter(f)),
                                                f += _,
                                                l += _,
                                                y -= 1;
                                        return r
                                    }
                                }),
                                Object.defineProperty(l.prototype, "every", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        for (var i = arguments[1], o = 0; o < n; o++)
                                            if (!t.call(i, e._getter(o), o, e))
                                                return !1;
                                        return !0
                                    }
                                }),
                                Object.defineProperty(l.prototype, "fill", {
                                    value: function(t) {
                                        var e = arguments[1]
                                            , n = arguments[2]
                                            , r = i(this)
                                            , s = r.length
                                            , u = a(s);
                                        u = I(u, 0);
                                        var l, c = o(e);
                                        l = c < 0 ? I(u + c, 0) : L(c, u);
                                        var f;
                                        f = n === P ? u : o(n);
                                        var h;
                                        for (h = f < 0 ? I(u + f, 0) : L(f, u); l < h; )
                                            r._setter(l, t),
                                                l += 1;
                                        return r
                                    }
                                }),
                                Object.defineProperty(l.prototype, "filter", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        for (var i = [], o = arguments[1], s = 0; s < n; s++) {
                                            var u = e._getter(s);
                                            t.call(o, u, s, e) && i.push(u)
                                        }
                                        return new this.constructor(i)
                                    }
                                }),
                                Object.defineProperty(l.prototype, "find", {
                                    value: function(t) {
                                        var e = i(this)
                                            , n = e.length
                                            , o = a(n);
                                        if (!r(t))
                                            throw TypeError();
                                        for (var s = arguments.length > 1 ? arguments[1] : P, u = 0; u < o; ) {
                                            var l = e._getter(u)
                                                , c = t.call(s, l, u, e);
                                            if (Boolean(c))
                                                return l;
                                            ++u
                                        }
                                        return P
                                    }
                                }),
                                Object.defineProperty(l.prototype, "findIndex", {
                                    value: function(t) {
                                        var e = i(this)
                                            , n = e.length
                                            , o = a(n);
                                        if (!r(t))
                                            throw TypeError();
                                        for (var s = arguments.length > 1 ? arguments[1] : P, u = 0; u < o; ) {
                                            var l = e._getter(u)
                                                , c = t.call(s, l, u, e);
                                            if (Boolean(c))
                                                return u;
                                            ++u
                                        }
                                        return -1
                                    }
                                }),
                                Object.defineProperty(l.prototype, "forEach", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        for (var i = arguments[1], o = 0; o < n; o++)
                                            t.call(i, e._getter(o), o, e)
                                    }
                                }),
                                Object.defineProperty(l.prototype, "indexOf", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (0 === n)
                                            return -1;
                                        var r = 0;
                                        if (arguments.length > 0 && (r = Number(arguments[1]),
                                            r !== r ? r = 0 : 0 !== r && r !== 1 / 0 && r !== -(1 / 0) && (r = (r > 0 || -1) * N(z(r)))),
                                        r >= n)
                                            return -1;
                                        for (var i = r >= 0 ? r : I(n - z(r), 0); i < n; i++)
                                            if (e._getter(i) === t)
                                                return i;
                                        return -1
                                    }
                                }),
                                Object.defineProperty(l.prototype, "join", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        for (var e = Object(this), n = a(e.length), r = Array(n), i = 0; i < n; ++i)
                                            r[i] = e._getter(i);
                                        return r.join(t === P ? "," : t)
                                    }
                                }),
                                Object.defineProperty(l.prototype, "lastIndexOf", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (0 === n)
                                            return -1;
                                        var r = n;
                                        arguments.length > 1 && (r = Number(arguments[1]),
                                            r !== r ? r = 0 : 0 !== r && r !== 1 / 0 && r !== -(1 / 0) && (r = (r > 0 || -1) * N(z(r))));
                                        for (var i = r >= 0 ? L(r, n - 1) : n - z(r); i >= 0; i--)
                                            if (e._getter(i) === t)
                                                return i;
                                        return -1
                                    }
                                }),
                                Object.defineProperty(l.prototype, "map", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        var i = [];
                                        i.length = n;
                                        for (var o = arguments[1], s = 0; s < n; s++)
                                            i[s] = t.call(o, e._getter(s), s, e);
                                        return new this.constructor(i)
                                    }
                                }),
                                Object.defineProperty(l.prototype, "reduce", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        if (0 === n && 1 === arguments.length)
                                            throw TypeError();
                                        var i, o = 0;
                                        for (i = arguments.length >= 2 ? arguments[1] : e._getter(o++); o < n; )
                                            i = t.call(P, i, e._getter(o), o, e),
                                                o++;
                                        return i
                                    }
                                }),
                                Object.defineProperty(l.prototype, "reduceRight", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        if (0 === n && 1 === arguments.length)
                                            throw TypeError();
                                        var i, o = n - 1;
                                        for (i = arguments.length >= 2 ? arguments[1] : e._getter(o--); o >= 0; )
                                            i = t.call(P, i, e._getter(o), o, e),
                                                o--;
                                        return i
                                    }
                                }),
                                Object.defineProperty(l.prototype, "reverse", {
                                    value: function() {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        for (var t = Object(this), e = a(t.length), n = N(e / 2), r = 0, i = e - 1; r < n; ++r,
                                            --i) {
                                            var o = t._getter(r);
                                            t._setter(r, t._getter(i)),
                                                t._setter(i, o)
                                        }
                                        return t
                                    }
                                }),
                                Object.defineProperty(l.prototype, "set", {
                                    value: function(t, e) {
                                        if (arguments.length < 1)
                                            throw SyntaxError("Not enough arguments");
                                        var n, r, i, o, s, u, l, c, f, h;
                                        if ("object" === babelHelpers.typeof(arguments[0]) && arguments[0].constructor === this.constructor) {
                                            if (n = arguments[0],
                                                i = a(arguments[1]),
                                            i + n.length > this.length)
                                                throw RangeError("Offset plus length of array is out of range");
                                            if (c = this.byteOffset + i * this.BYTES_PER_ELEMENT,
                                                f = n.length * this.BYTES_PER_ELEMENT,
                                            n.buffer === this.buffer) {
                                                for (h = [],
                                                        s = 0,
                                                        u = n.byteOffset; s < f; s += 1,
                                                        u += 1)
                                                    h[s] = n.buffer._bytes[u];
                                                for (s = 0,
                                                        l = c; s < f; s += 1,
                                                        l += 1)
                                                    this.buffer._bytes[l] = h[s]
                                            } else
                                                for (s = 0,
                                                        u = n.byteOffset,
                                                        l = c; s < f; s += 1,
                                                        u += 1,
                                                        l += 1)
                                                    this.buffer._bytes[l] = n.buffer._bytes[u]
                                        } else {
                                            if ("object" !== babelHelpers.typeof(arguments[0]) || "undefined" == typeof arguments[0].length)
                                                throw TypeError("Unexpected argument type(s)");
                                            if (r = arguments[0],
                                                o = a(r.length),
                                                i = a(arguments[1]),
                                            i + o > this.length)
                                                throw RangeError("Offset plus length of array is out of range");
                                            for (s = 0; s < o; s += 1)
                                                u = r[s],
                                                    this._setter(i + s, Number(u))
                                        }
                                    }
                                }),
                                Object.defineProperty(l.prototype, "slice", {
                                    value: function(t, e) {
                                        for (var n = i(this), r = n.length, s = a(r), u = o(t), l = u < 0 ? I(s + u, 0) : L(u, s), c = e === P ? s : o(e), f = c < 0 ? I(s + c, 0) : L(c, s), h = f - l, d = n.constructor, p = new d(h), _ = 0; l < f; ) {
                                            var y = n._getter(l);
                                            p._setter(_, y),
                                                ++l,
                                                ++_
                                        }
                                        return p
                                    }
                                }),
                                Object.defineProperty(l.prototype, "some", {
                                    value: function(t) {
                                        if (this === P || null === this)
                                            throw TypeError();
                                        var e = Object(this)
                                            , n = a(e.length);
                                        if (!r(t))
                                            throw TypeError();
                                        for (var i = arguments[1], o = 0; o < n; o++)
                                            if (t.call(i, e._getter(o), o, e))
                                                return !0;
                                        return !1
                                    }
                                }),
                                Object.defineProperty(l.prototype, "sort", {
                                    value: function(t) {
                                        function e(e, n) {
                                            return e !== e && n !== n ? 0 : e !== e ? 1 : n !== n ? -1 : t !== P ? t(e, n) : e < n ? -1 : e > n ? 1 : 0
                                        }
                                        if (this === P || null === this)
                                            throw TypeError();
                                        for (var n = Object(this), r = a(n.length), i = Array(r), o = 0; o < r; ++o)
                                            i[o] = n._getter(o);
                                        for (i.sort(e),
                                                o = 0; o < r; ++o)
                                            n._setter(o, i[o]);
                                        return n
                                    }
                                }),
                                Object.defineProperty(l.prototype, "subarray", {
                                    value: function(t, e) {
                                        function n(t, e, n) {
                                            return t < e ? e : t > n ? n : t
                                        }
                                        t = o(t),
                                            e = o(e),
                                        arguments.length < 1 && (t = 0),
                                        arguments.length < 2 && (e = this.length),
                                        t < 0 && (t = this.length + t),
                                        e < 0 && (e = this.length + e),
                                            t = n(t, 0, this.length),
                                            e = n(e, 0, this.length);
                                        var r = e - t;
                                        return r < 0 && (r = 0),
                                            new this.constructor(this.buffer,this.byteOffset + t * this.BYTES_PER_ELEMENT,r)
                                    }
                                });
                            var C = k(1, c, f)
                                , S = k(1, h, d)
                                , R = k(1, p, d)
                                , B = k(2, _, y)
                                , M = k(2, b, g)
                                , F = k(4, v, w)
                                , Y = k(4, m, E)
                                , U = k(4, x, A)
                                , D = k(8, j, T)
                                , X = document.documentMode || +(navigator.userAgent.match(/MSIE (\d+)/) && RegExp.$1) || !t.Int8Array;
                            t.Int8Array = X ? C : t.Int8Array,
                                t.Uint8Array = X ? S : t.Uint8Array,
                                t.Uint8ClampedArray = X ? R : t.Uint8ClampedArray,
                                t.Int16Array = X ? B : t.Int16Array,
                                t.Uint16Array = X ? M : t.Uint16Array,
                                t.Int32Array = X ? F : t.Int32Array,
                                t.Uint32Array = X ? Y : t.Uint32Array,
                                t.Float32Array = X ? U : t.Float32Array,
                                t.Float64Array = X ? D : t.Float64Array
                        }(),
                        function() {
                            function e(t, e) {
                                return r(t.get) ? t.get(e) : t[e]
                            }
                            function i(t, e, r) {
                                if (!(t instanceof ArrayBuffer || "ArrayBuffer" === n(t)))
                                    throw TypeError();
                                if (e = a(e),
                                e > t.byteLength)
                                    throw RangeError("byteOffset out of range");
                                if (r = r === P ? t.byteLength - e : a(r),
                                e + r > t.byteLength)
                                    throw RangeError("byteOffset and length reference an area beyond the end of the buffer");
                                Object.defineProperty(this, "buffer", {
                                    value: t
                                }),
                                    Object.defineProperty(this, "byteLength", {
                                        value: r
                                    }),
                                    Object.defineProperty(this, "byteOffset", {
                                        value: e
                                    })
                            }
                            function o(t) {
                                return function(n, r) {
                                    if (n = a(n),
                                    n + t.BYTES_PER_ELEMENT > this.byteLength)
                                        throw RangeError("Array index out of range");
                                    n += this.byteOffset;
                                    for (var i = new Uint8Array(this.buffer,n,t.BYTES_PER_ELEMENT), o = [], s = 0; s < t.BYTES_PER_ELEMENT; s += 1)
                                        o.push(e(i, s));
                                    return Boolean(r) === Boolean(u) && o.reverse(),
                                        e(new t(new Uint8Array(o).buffer), 0)
                                }
                            }
                            function s(t) {
                                return function(n, r, i) {
                                    if (n = a(n),
                                    n + t.BYTES_PER_ELEMENT > this.byteLength)
                                        throw RangeError("Array index out of range");
                                    var o, s, l = new t([r]), c = new Uint8Array(l.buffer), f = [];
                                    for (o = 0; o < t.BYTES_PER_ELEMENT; o += 1)
                                        f.push(e(c, o));
                                    Boolean(i) === Boolean(u) && f.reverse(),
                                        s = new Uint8Array(this.buffer,n,t.BYTES_PER_ELEMENT),
                                        s.set(f)
                                }
                            }
                            var u = function() {
                                var t = new Uint16Array([4660])
                                    , n = new Uint8Array(t.buffer);
                                return 18 === e(n, 0)
                            }();
                            Object.defineProperty(i.prototype, "getUint8", {
                                value: o(Uint8Array)
                            }),
                                Object.defineProperty(i.prototype, "getInt8", {
                                    value: o(Int8Array)
                                }),
                                Object.defineProperty(i.prototype, "getUint16", {
                                    value: o(Uint16Array)
                                }),
                                Object.defineProperty(i.prototype, "getInt16", {
                                    value: o(Int16Array)
                                }),
                                Object.defineProperty(i.prototype, "getUint32", {
                                    value: o(Uint32Array)
                                }),
                                Object.defineProperty(i.prototype, "getInt32", {
                                    value: o(Int32Array)
                                }),
                                Object.defineProperty(i.prototype, "getFloat32", {
                                    value: o(Float32Array)
                                }),
                                Object.defineProperty(i.prototype, "getFloat64", {
                                    value: o(Float64Array)
                                }),
                                Object.defineProperty(i.prototype, "setUint8", {
                                    value: s(Uint8Array)
                                }),
                                Object.defineProperty(i.prototype, "setInt8", {
                                    value: s(Int8Array)
                                }),
                                Object.defineProperty(i.prototype, "setUint16", {
                                    value: s(Uint16Array)
                                }),
                                Object.defineProperty(i.prototype, "setInt16", {
                                    value: s(Int16Array)
                                }),
                                Object.defineProperty(i.prototype, "setUint32", {
                                    value: s(Uint32Array)
                                }),
                                Object.defineProperty(i.prototype, "setInt32", {
                                    value: s(Int32Array)
                                }),
                                Object.defineProperty(i.prototype, "setFloat32", {
                                    value: s(Float32Array)
                                }),
                                Object.defineProperty(i.prototype, "setFloat64", {
                                    value: s(Float64Array)
                                }),
                                t.DataView = t.DataView || i
                        }()
                }
        }
        , function(t, e, n) {
            "use strict";
            function r(t, e, n, r) {
                return r = r || {},
                    r["Content-Type"] = "application/x-www-form-urlencoded",
                    new u.default(function(i, o) {
                            var a = Date.now()
                                , s = new XMLHttpRequest;
                            if ("withCredentials"in s) {
                                if (s.open(e, t, !0),
                                    r)
                                    for (var u in r)
                                        r.hasOwnProperty(u) && s.setRequestHeader(u, r[u]);
                                s.onload = function() {
                                    if (4 === s.readyState)
                                        if (s.status >= 200 && s.status < 300) {
                                            var e = JSON.parse(s.response);
                                            window.Yoda && window.Yoda.CAT && window.Yoda.CAT.postBatch(t, 0, 0, Date.now() - a, "200|" + e.status, "ajax"),
                                                i(e),
                                                s = null
                                        } else
                                            o(new Error(s.statusText)),
                                                s = null
                                }
                                    ,
                                    s.ontimeout = function(e) {
                                        window.Yoda.CAT.postBatch(t, 0, 0, Date.now() - a, "500|0", "ajax"),
                                            window.Yoda.CAT.sendLog(t, "ajaxError", "【请求后端API接口超时】", e.message),
                                            o(new Error("请求超时:" + t)),
                                            s = null
                                    }
                                    ,
                                    s.onerror = function(e) {
                                        window.Yoda.CAT.postBatch(t, 0, 0, Date.now() - a, "500|0", "ajax"),
                                            window.Yoda.CAT.sendLog(t, "ajaxError", "【请求后端API接口ERROR】", e.message),
                                            o(new Error(s.statusText)),
                                            s = null
                                    }
                                    ,
                                    s.send(n)
                            } else
                                "undefined" != typeof XDomainRequest ? (0,
                                    c.default)({
                                    url: t,
                                    callback: "callback",
                                    data: n,
                                    success: function(t) {
                                        i(t)
                                    },
                                    fail: function(t) {
                                        o(new Error(t))
                                    },
                                    time: 1e4
                                }) : (o(new Error("创建xhr对象失败")),
                                    window.Yoda.CAT.sendLog(t, "ajaxError", "【请求后端API接口】", "创建XMLHttpRequest对象失败"),
                                    s = null)
                        }
                    ).catch(function(e) {
                        "production" === window.YODA_CONFIG.__ENV__ && window.Yoda.CAT.sendLog(t, "ajaxError", "【请求后端API接口】:发生异常promise-catch", e.message)
                    })
            }
            function i(t) {
                var e = "&";
                return t.indexOf("?") === -1 && (e = "?"),
                    t += e + a("GET", t, ""),
                    r(t, "GET", null)
            }
            function o(t, e, n) {
                if (null !== e && "object" === ("undefined" == typeof e ? "undefined" : babelHelpers.typeof(e)) && !(e instanceof String || window.FormData && e instanceof window.FormData)) {
                    var i = [];
                    for (var o in e)
                        i.push(encodeURIComponent(o) + "=" + encodeURIComponent(e[o]));
                    e = i.join("&")
                }
                var s = "&";
                return (!e || e.length < 1) && (s = ""),
                    e += s + a("POST", t, e),
                    r(t, "POST", e, n)
            }
            function a(t, e, n) {
                try {
                    if (e.indexOf(d._a) > 0 || n.length > 0 && n.indexOf(d._a) > 0)
                        return "";
                    var r = "";
                    return "GET" === t ? r = d.reload(e) : (e = e.indexOf("?") > 0 ? e + "&" + n : e + "?" + n,
                        r = d.reload(e)),
                    r || window.Yoda.CAT.sendLog(e, "jsError", "【url参数处理异常】", "t 为空"),
                    encodeURIComponent(d._a) + "=" + encodeURIComponent(r)
                } catch (t) {
                    "production" === window.YODA_CONFIG.__ENV__ && window.Yoda.CAT.sendLog(e, "jsError", "【url参数处理异常】", t.message)
                }
            }
            var s = n(6)
                , u = babelHelpers.interopRequireDefault(s)
                , l = n(10)
                , c = babelHelpers.interopRequireDefault(l)
                , f = n(11)
                , h = babelHelpers.interopRequireDefault(f)
                , d = new h.default
                , p = {
                post: o,
                get: i
            };
            t.exports = p
        }
        , function(t, e, n) {
            (function(e) {
                    !function(n) {
                        function r() {}
                        function i(t, e) {
                            return function() {
                                t.apply(e, arguments)
                            }
                        }
                        function o(t) {
                            if ("object" != typeof this)
                                throw new TypeError("Promises must be constructed via new");
                            if ("function" != typeof t)
                                throw new TypeError("not a function");
                            this._state = 0,
                                this._handled = !1,
                                this._value = void 0,
                                this._deferreds = [],
                                f(t, this)
                        }
                        function a(t, e) {
                            for (; 3 === t._state; )
                                t = t._value;
                            return 0 === t._state ? void t._deferreds.push(e) : (t._handled = !0,
                                void o._immediateFn(function() {
                                    var n = 1 === t._state ? e.onFulfilled : e.onRejected;
                                    if (null === n)
                                        return void (1 === t._state ? s : u)(e.promise, t._value);
                                    var r;
                                    try {
                                        r = n(t._value)
                                    } catch (t) {
                                        return void u(e.promise, t)
                                    }
                                    s(e.promise, r)
                                }))
                        }
                        function s(t, e) {
                            try {
                                if (e === t)
                                    throw new TypeError("A promise cannot be resolved with itself.");
                                if (e && ("object" == typeof e || "function" == typeof e)) {
                                    var n = e.then;
                                    if (e instanceof o)
                                        return t._state = 3,
                                            t._value = e,
                                            void l(t);
                                    if ("function" == typeof n)
                                        return void f(i(n, e), t)
                                }
                                t._state = 1,
                                    t._value = e,
                                    l(t)
                            } catch (e) {
                                u(t, e)
                            }
                        }
                        function u(t, e) {
                            t._state = 2,
                                t._value = e,
                                l(t)
                        }
                        function l(t) {
                            2 === t._state && 0 === t._deferreds.length && o._immediateFn(function() {
                                t._handled || o._unhandledRejectionFn(t._value)
                            });
                            for (var e = 0, n = t._deferreds.length; e < n; e++)
                                a(t, t._deferreds[e]);
                            t._deferreds = null
                        }
                        function c(t, e, n) {
                            this.onFulfilled = "function" == typeof t ? t : null,
                                this.onRejected = "function" == typeof e ? e : null,
                                this.promise = n
                        }
                        function f(t, e) {
                            var n = !1;
                            try {
                                t(function(t) {
                                    n || (n = !0,
                                        s(e, t))
                                }, function(t) {
                                    n || (n = !0,
                                        u(e, t))
                                })
                            } catch (t) {
                                if (n)
                                    return;
                                n = !0,
                                    u(e, t)
                            }
                        }
                        var h = setTimeout;
                        o.prototype.catch = function(t) {
                            return this.then(null, t)
                        }
                            ,
                            o.prototype.then = function(t, e) {
                                var n = new this.constructor(r);
                                return a(this, new c(t,e,n)),
                                    n
                            }
                            ,
                            o.all = function(t) {
                                var e = Array.prototype.slice.call(t);
                                return new o(function(t, n) {
                                        function r(o, a) {
                                            try {
                                                if (a && ("object" == typeof a || "function" == typeof a)) {
                                                    var s = a.then;
                                                    if ("function" == typeof s)
                                                        return void s.call(a, function(t) {
                                                            r(o, t)
                                                        }, n)
                                                }
                                                e[o] = a,
                                                0 === --i && t(e)
                                            } catch (t) {
                                                n(t)
                                            }
                                        }
                                        if (0 === e.length)
                                            return t([]);
                                        for (var i = e.length, o = 0; o < e.length; o++)
                                            r(o, e[o])
                                    }
                                )
                            }
                            ,
                            o.resolve = function(t) {
                                return t && "object" == typeof t && t.constructor === o ? t : new o(function(e) {
                                        e(t)
                                    }
                                )
                            }
                            ,
                            o.reject = function(t) {
                                return new o(function(e, n) {
                                        n(t)
                                    }
                                )
                            }
                            ,
                            o.race = function(t) {
                                return new o(function(e, n) {
                                        for (var r = 0, i = t.length; r < i; r++)
                                            t[r].then(e, n)
                                    }
                                )
                            }
                            ,
                            o._immediateFn = "function" == typeof e && function(t) {
                                    e(t)
                                }
                                || function(t) {
                                    h(t, 0)
                                }
                            ,
                            o._unhandledRejectionFn = function(t) {
                                "undefined" != typeof console && console
                            }
                            ,
                            o._setImmediateFn = function(t) {
                                o._immediateFn = t
                            }
                            ,
                            o._setUnhandledRejectionFn = function(t) {
                                o._unhandledRejectionFn = t
                            }
                            ,
                            "undefined" != typeof t && t.exports ? t.exports = o : n.Promise || (n.Promise = o)
                    }(this)
                }
            ).call(e, n(7).setImmediate)
        }
        , function(t, e, n) {
            function r(t, e) {
                this._id = t,
                    this._clearFn = e
            }
            var i = Function.prototype.apply;
            e.setTimeout = function() {
                return new r(i.call(setTimeout, window, arguments),clearTimeout)
            }
                ,
                e.setInterval = function() {
                    return new r(i.call(setInterval, window, arguments),clearInterval)
                }
                ,
                e.clearTimeout = e.clearInterval = function(t) {
                    t && t.close()
                }
                ,
                r.prototype.unref = r.prototype.ref = function() {}
                ,
                r.prototype.close = function() {
                    this._clearFn.call(window, this._id)
                }
                ,
                e.enroll = function(t, e) {
                    clearTimeout(t._idleTimeoutId),
                        t._idleTimeout = e
                }
                ,
                e.unenroll = function(t) {
                    clearTimeout(t._idleTimeoutId),
                        t._idleTimeout = -1
                }
                ,
                e._unrefActive = e.active = function(t) {
                    clearTimeout(t._idleTimeoutId);
                    var e = t._idleTimeout;
                    e >= 0 && (t._idleTimeoutId = setTimeout(function() {
                        t._onTimeout && t._onTimeout()
                    }, e))
                }
                ,
                n(8),
                e.setImmediate = setImmediate,
                e.clearImmediate = clearImmediate
        }
        , function(t, e, n) {
            (function(t, e) {
                    !function(t, n) {
                        "use strict";
                        function r(t) {
                            "function" != typeof t && (t = new Function("" + t));
                            for (var e = new Array(arguments.length - 1), n = 0; n < e.length; n++)
                                e[n] = arguments[n + 1];
                            var r = {
                                callback: t,
                                args: e
                            };
                            return _[p] = r,
                                d(p),
                                p++
                        }
                        function i(t) {
                            delete _[t]
                        }
                        function o(t) {
                            var e = t.callback
                                , r = t.args;
                            switch (r.length) {
                                case 0:
                                    e();
                                    break;
                                case 1:
                                    e(r[0]);
                                    break;
                                case 2:
                                    e(r[0], r[1]);
                                    break;
                                case 3:
                                    e(r[0], r[1], r[2]);
                                    break;
                                default:
                                    e.apply(n, r)
                            }
                        }
                        function a(t) {
                            if (y)
                                setTimeout(a, 0, t);
                            else {
                                var e = _[t];
                                if (e) {
                                    y = !0;
                                    try {
                                        o(e)
                                    } finally {
                                        i(t),
                                            y = !1
                                    }
                                }
                            }
                        }
                        function s() {
                            d = function(t) {
                                e.nextTick(function() {
                                    a(t)
                                })
                            }
                        }
                        function u() {
                            if (t.postMessage && !t.importScripts) {
                                var e = !0
                                    , n = t.onmessage;
                                return t.onmessage = function() {
                                    e = !1
                                }
                                    ,
                                    t.postMessage("", "*"),
                                    t.onmessage = n,
                                    e
                            }
                        }
                        function l() {
                            var e = "setImmediate$" + Math.random() + "$"
                                , n = function(n) {
                                n.source === t && "string" == typeof n.data && 0 === n.data.indexOf(e) && a(+n.data.slice(e.length))
                            };
                            t.addEventListener ? t.addEventListener("message", n, !1) : t.attachEvent("onmessage", n),
                                d = function(n) {
                                    t.postMessage(e + n, "*")
                                }
                        }
                        function c() {
                            var t = new MessageChannel;
                            t.port1.onmessage = function(t) {
                                var e = t.data;
                                a(e)
                            }
                                ,
                                d = function(e) {
                                    t.port2.postMessage(e)
                                }
                        }
                        function f() {
                            var t = b.documentElement;
                            d = function(e) {
                                var n = b.createElement("script");
                                n.onreadystatechange = function() {
                                    a(e),
                                        n.onreadystatechange = null,
                                        t.removeChild(n),
                                        n = null
                                }
                                    ,
                                    t.appendChild(n)
                            }
                        }
                        function h() {
                            d = function(t) {
                                setTimeout(a, 0, t)
                            }
                        }
                        if (!t.setImmediate) {
                            var d, p = 1, _ = {}, y = !1, b = t.document, g = Object.getPrototypeOf && Object.getPrototypeOf(t);
                            g = g && g.setTimeout ? g : t,
                                "[object process]" === {}.toString.call(t.process) ? s() : u() ? l() : t.MessageChannel ? c() : b && "onreadystatechange"in b.createElement("script") ? f() : h(),
                                g.setImmediate = r,
                                g.clearImmediate = i
                        }
                    }("undefined" == typeof self ? "undefined" == typeof t ? this : t : self)
                }
            ).call(e, function() {
                return this
            }(), n(9))
        }
        , function(t, e) {
            function n() {
                throw new Error("setTimeout has not been defined")
            }
            function r() {
                throw new Error("clearTimeout has not been defined")
            }
            function i(t) {
                if (c === setTimeout)
                    return setTimeout(t, 0);
                if ((c === n || !c) && setTimeout)
                    return c = setTimeout,
                        setTimeout(t, 0);
                try {
                    return c(t, 0)
                } catch (e) {
                    try {
                        return c.call(null, t, 0)
                    } catch (e) {
                        return c.call(this, t, 0)
                    }
                }
            }
            function o(t) {
                if (f === clearTimeout)
                    return clearTimeout(t);
                if ((f === r || !f) && clearTimeout)
                    return f = clearTimeout,
                        clearTimeout(t);
                try {
                    return f(t)
                } catch (e) {
                    try {
                        return f.call(null, t)
                    } catch (e) {
                        return f.call(this, t)
                    }
                }
            }
            function a() {
                _ && d && (_ = !1,
                    d.length ? p = d.concat(p) : y = -1,
                p.length && s())
            }
            function s() {
                if (!_) {
                    var t = i(a);
                    _ = !0;
                    for (var e = p.length; e; ) {
                        for (d = p,
                                p = []; ++y < e; )
                            d && d[y].run();
                        y = -1,
                            e = p.length
                    }
                    d = null,
                        _ = !1,
                        o(t)
                }
            }
            function u(t, e) {
                this.fun = t,
                    this.array = e
            }
            function l() {}
            var c, f, h = t.exports = {};
            !function() {
                try {
                    c = "function" == typeof setTimeout ? setTimeout : n
                } catch (t) {
                    c = n
                }
                try {
                    f = "function" == typeof clearTimeout ? clearTimeout : r
                } catch (t) {
                    f = r
                }
            }();
            var d, p = [], _ = !1, y = -1;
            h.nextTick = function(t) {
                var e = new Array(arguments.length - 1);
                if (arguments.length > 1)
                    for (var n = 1; n < arguments.length; n++)
                        e[n - 1] = arguments[n];
                p.push(new u(t,e)),
                1 !== p.length || _ || i(s)
            }
                ,
                u.prototype.run = function() {
                    this.fun.apply(null, this.array)
                }
                ,
                h.title = "browser",
                h.browser = !0,
                h.env = {},
                h.argv = [],
                h.version = "",
                h.versions = {},
                h.on = l,
                h.addListener = l,
                h.once = l,
                h.off = l,
                h.removeListener = l,
                h.removeAllListeners = l,
                h.emit = l,
                h.prependListener = l,
                h.prependOnceListener = l,
                h.listeners = function(t) {
                    return []
                }
                ,
                h.binding = function(t) {
                    throw new Error("process.binding is not supported")
                }
                ,
                h.cwd = function() {
                    return "/"
                }
                ,
                h.chdir = function(t) {
                    throw new Error("process.chdir is not supported")
                }
                ,
                h.umask = function() {
                    return 0
                }
        }
        , function(t, e) {
            "use strict";
            Object.defineProperty(e, "__esModule", {
                value: !0
            });
            var n = function(t) {
                if (t = t || {},
                !t.url || !t.callback)
                    throw new Error("参数异常,请检查");
                var e = ("jsonp_" + Math.random()).replace(".", "")
                    , n = document.getElementsByTagName("head")[0]
                    , r = "";
                t.data ? r = t.data + "&" + t.callback + "=" + e : r += t.callback + "=" + e;
                var i = document.createElement("script");
                n.appendChild(i),
                    window[e] = function(r) {
                        n.removeChild(i),
                            clearTimeout(i.timer),
                            window[e] = null,
                        t.success && t.success(r)
                    }
                    ,
                    t.url.indexOf("?") ? i.src = t.url + "&" + r : i.src = t.url + "?" + r,
                t.time && (i.timer = setTimeout(function() {
                    window[e] = null,
                        n.removeChild(i),
                    t.fail && t.fail({
                        message: "请求超时"
                    })
                }, t.time))
            };
            e.default = n
        }
        , function(t, e, n) {
            "use strict";
            var r = function t() {
                var e = n(12).deflate
                    , r = n(21)
                    , i = n(24);
                Object.keys || (Object.keys = n(25)),
                Function.prototype.bind || (Function.prototype.bind = function(t) {
                        if ("function" != typeof this)
                            throw new TypeError("Function.prototype.bind - what is trying to be bound is not callable");
                        var e = Array.prototype.slice.call(arguments, 1)
                            , n = this
                            , r = function() {}
                            , i = function() {
                            return n.apply(this instanceof r && t ? this : t, e.concat(Array.prototype.slice.call(arguments)))
                        };
                        return r.prototype = this.prototype,
                            i.prototype = new r,
                            i
                    }
                ),
                "function" != typeof Array.prototype.forEach && (Array.prototype.forEach = function(t, e) {
                        for (var n = 0; n < this.length; n++)
                            t.apply(e, [this[n], n, this])
                    }
                ),
                "undefined" == typeof JSON && (JSON = n(27));
                var o = function() {
                    var t = Math.max(document.documentElement.clientWidth, window.innerWidth || 0)
                        , e = Math.max(document.documentElement.clientHeight, window.innerHeight || 0);
                    return [t, e]
                }
                    , a = function() {
                    var t = [screen.width, screen.height]
                        , e = [screen.availWidth, screen.availHeight]
                        , n = screen.colorDepth
                        , r = screen.pixelDepth;
                    return [t, e, n, r]
                }
                    , s = function() {
                    try {
                        var t = Function("return this")()
                            , e = function() {
                            var e = (t.constructor + "").match(/ (\w+)|$/)[1];
                            if (!e)
                                try {
                                    "[object]" == t && (e = "Window")
                                } catch (t) {
                                    e = "WSH"
                                }
                            return e
                        }()
                            , n = "";
                        switch (e) {
                            case "Window":
                                break;
                            case "DedicatedWorkerGlobalScope":
                                n = "ww";
                                break;
                            case "WSH":
                                n = "wsh";
                                break;
                            case "Object":
                                n = "nj";
                                break;
                            default:
                                n = "ot"
                        }
                        return n
                    } catch (t) {
                        return "abnormal"
                    }
                }
                    , u = function() {
                    return window._phantom || window.phantom || window.callPhantom ? "ps" : s() || i.getwd()
                }
                    , l = function() {
                    var t = document.referrer
                        , e = window.location.href;
                    return [e, t]
                }
                    , c = function(t) {
                    try {
                        t = e(JSON.stringify(t), {
                            to: "string"
                        })
                    } catch (e) {
                        throw new Error(t + " - 错误信息:" + e.message)
                    }
                    try {
                        t = btoa(t)
                    } catch (e) {
                        throw new Error(t + " - 错误信息:" + e.message)
                    }
                    return t
                }
                    , f = function(e) {
                    var n = []
                        , r = Object.keys(e).sort();
                    return r.forEach(function(r, i) {
                        r !== t._a && n.push(r + "=" + e[r])
                    }),
                        n = n.join("&"),
                        c(n)
                }
                    , h = function(t) {
                    t = t || window.event;
                    var e = t.pageX || t.clientX + (document.documentElement.scrollLeft || document.body.scrollLeft)
                        , n = t.pageY || t.clientY + (document.documentElement.scrollTop || document.body.scrollTop);
                    return {
                        x: e,
                        y: n
                    }
                }
                    , d = function() {
                    var t, e = window.navigator, n = e.plugins, r = [];
                    for (t in n)
                        if (n.hasOwnProperty(t)) {
                            var i = n[t].name || "";
                            r.push(i)
                        }
                    return r
                }
                    , t = {
                    v: "2.0.1",
                    rId: "100009",
                    ts: (new Date).getTime(),
                    cts: (new Date).getTime(),
                    brVD: o(),
                    brR: a(),
                    bI: l(),
                    mT: [],
                    kT: [],
                    aT: [],
                    tT: [],
                    aM: u(),
                    inputs: [],
                    buttons: [],
                    broP: d()
                };
                return t.bindUserTrackEvent = function() {
                    function e(t, e, n, r) {
                        e.addEventListener ? e.addEventListener(t, n, r || !1) : e.attachEvent ? e.attachEvent("on" + t, n) : e[t] = n
                    }
                    var n = function(e) {
                        var n, r, i;
                        e = e || window.event,
                        null == e.pageX && null !== e.clientX && (n = e.target && e.target.ownerDocument || document,
                            r = n.documentElement,
                            i = n.body,
                            e.pageX = e.clientX + (r && r.scrollLeft || i && i.scrollLeft || 0) - (r && r.clientLeft || i && i.clientLeft || 0),
                            e.pageY = e.clientY + (r && r.scrollTop || i && i.scrollTop || 0) - (r && r.clientTop || i && i.clientTop || 0));
                        var o = Date.now() - t.ts;
                        this.mT.unshift([e.pageX, e.pageY, o].join(",")),
                        this.mT.length > 60 && (this.mT = this.mT.slice(0, 60))
                    }
                        .bind(this)
                        , r = function(e) {
                        e = e || window.event;
                        var n = e.target || e.srcElement
                            , r = "number" == typeof e.which ? e.which : e.keyCode;
                        if (r) {
                            var i = Date.now() - t.ts;
                            this.kT.unshift([String.fromCharCode(r), n.nodeName, i].join(","))
                        }
                        this.kT.length > 30 && (this.kT = this.kT.slice(0, 30))
                    }
                        .bind(this)
                        , o = function(e) {
                        var n, r, i, o, a;
                        null !== e.touches[0].clientX && (n = e.target && e.target.ownerDocument || document,
                            r = n.documentElement,
                            i = n.body,
                            o = e.touches[0].clientX + (r && r.scrollLeft || i && i.scrollLeft || 0) - (r && r.clientLeft || i && i.clientLeft || 0),
                            a = e.touches[0].clientY + (r && r.scrollTop || i && i.scrollTop || 0) - (r && r.clientTop || i && i.clientTop || 0));
                        var s = Date.now() - t.ts;
                        this.tT.unshift([o, a, e.touches.length, s].join(",")),
                        this.tT.length > 60 && (this.tT = this.tT.slice(0, 60))
                    }
                        .bind(this)
                        , a = function(e) {
                        e = e || window.event;
                        var n = e.target || e.srcElement
                            , r = Date.now() - t.ts;
                        this.aT.unshift([e.clientX, e.clientY, n.nodeName, r].join(",")),
                        this.aT.length > 30 && (this.aT = this.aT.slice(0, 30))
                    }
                        .bind(this);
                    e("mousemove", document, n, !0),
                        e("keydown", document, r, !0),
                        e("click", document, a, !0),
                    "ontouchmove"in document && e("touchmove", document, o, !0),
                    0 === t.aM.length && i.listenwd(function(e) {
                        e && e.length > 0 && (t.aM = e)
                    });
                    var s = function(t) {
                        t = t || window.event;
                        var e = t.target || t.srcElement;
                        if (e.nodeName && "INPUT" === e.nodeName) {
                            var n = e.name || e.id;
                            n || (n = e.id = "rohr_" + parseInt(1e6 * Math.random()));
                            for (var r = this.inputs.length, i = 0; i < r; i++)
                                n === this.inputs[0].inputName && (this.inputs.splice(0, 1),
                                    i = 0,
                                    r -= 1);
                            this.inputs.unshift({
                                inputName: n,
                                editStartedTimeStamp: Date.now(),
                                keyboardEvent: "0-0-0-0"
                            })
                        }
                    }
                        .bind(this)
                        , u = function(t) {
                        t = t || window.event;
                        var e = t.target || t.srcElement;
                        if (e.nodeName && "INPUT" === e.nodeName) {
                            var n = this.inputs[0];
                            if (n) {
                                var r = n.keyboardEvent.split("-");
                                r[2] = 1,
                                    n.keyboardEvent = r.join("-")
                            }
                        }
                    }
                        .bind(this)
                        , l = function(t) {
                        t = t || window.event;
                        var e = t.target || t.srcElement;
                        if (e.nodeName && "INPUT" === e.nodeName) {
                            var n = this.inputs[0]
                                , r = n.keyboardEvent.split("-")
                                , i = "number" == typeof t.which ? t.which : t.keyCode;
                            9 === i && (r[0] = parseInt(r[0]) + 1),
                                r[1] = parseInt(r[1]) + 1;
                            var o = Date.now();
                            if (n.lastTime) {
                                var a = n.lastTime;
                                r[3] = r[3] + "|" + parseInt(o - a, 36)
                            }
                            this.inputs[0].lastTime = o,
                                this.inputs[0].keyboardEvent = r.join("-")
                        }
                    }
                        .bind(this)
                        , c = function(t) {
                        t = t || window.event;
                        var e = t.target || t.srcElement;
                        if (e.nodeName && "INPUT" === e.nodeName) {
                            var n = this.inputs[0];
                            if (!n)
                                return;
                            n.editFinishedTimeStamp = Date.now();
                            var r = n.keyboardEvent.split("-");
                            0 != r[3] && (r[3] = r[3].substr(2)),
                                delete n.lastTime,
                                n.keyboardEvent = r.join("-")
                        }
                    }
                        .bind(this)
                        , f = function(t) {
                        t = t || window.event;
                        var e = t.target || t.srcElement;
                        if (e.nodeName && "BUTTON" === e.nodeName) {
                            var n = e.name || e.id;
                            n || (n = e.id = "rohr_" + parseInt(1e6 * Math.random()));
                            for (var r = this.buttons.length, i = 0; i < r; i++)
                                n === this.buttons[i].buttonName && (this.buttons.splice(i, 1),
                                    i = 0,
                                    r -= 1);
                            var o = h(t)
                                , a = e.clientWidth
                                , s = e.clientHeight
                                , u = t.offsetX / a * 1e3
                                , l = (s - t.offsetY) / s * 1e3;
                            this.buttons.unshift({
                                buttonName: n,
                                touchPoint: "{" + o.x + "," + o.y + "}",
                                touchPosition: "{" + Math.floor(u) / 10 + "," + Math.floor(l) / 10 + "}",
                                touchTimeStamp: Date.now()
                            })
                        }
                    }
                        .bind(this);
                    e("focus", document, s, !0),
                        e("mouseout", document, u, !0),
                        e("keydown", document, l, !0),
                        e("blur", document, c, !0),
                        "ontouchstart"in document ? e("touchstart", document, f, !0) : e("click", document, f, !0)
                }
                    ,
                    t.reload = function(e) {
                        var n, i = {};
                        return "string" == typeof e ? i = r.parse(e.split("?")[1]) : null !== e && "object" === ("undefined" == typeof e ? "undefined" : babelHelpers.typeof(e)) && (i = e),
                            t.sign = f(i),
                            t.cts = Date.now(),
                            n = c(t)
                    }
                    ,
                    Object.defineProperty(t, "_a", {
                        get: function() {
                            var t = ""
                                , e = 0
                                , n = 3;
                            for (e; e < 6; ) {
                                var r = "";
                                switch (n) {
                                    case 47:
                                        r = "e",
                                            n = 78;
                                        break;
                                    case 3:
                                        r = "_",
                                            n = 9;
                                        break;
                                    case 78:
                                        r = "n";
                                        break;
                                    case 9:
                                        r = "t",
                                            n = 36;
                                        break;
                                    case 36:
                                        r = "o",
                                            n = 5;
                                        break;
                                    default:
                                        n = 47,
                                            r = "k"
                                }
                                e++,
                                    t += r
                            }
                            return t
                        }
                    }),
                    t.bindUserTrackEvent(),
                    {
                        reload: t.reload,
                        _a: t._a
                    }
            };
            t.exports = r
        }
        , function(t, e, n) {
            "use strict";
            function r(t) {
                if (!(this instanceof r))
                    return new r(t);
                this.options = u.assign({
                    level: g,
                    method: w,
                    chunkSize: 16384,
                    windowBits: 15,
                    memLevel: 8,
                    strategy: v,
                    to: ""
                }, t || {});
                var e = this.options;
                e.raw && e.windowBits > 0 ? e.windowBits = -e.windowBits : e.gzip && e.windowBits > 0 && e.windowBits < 16 && (e.windowBits += 16),
                    this.err = 0,
                    this.msg = "",
                    this.ended = !1,
                    this.chunks = [],
                    this.strm = new f,
                    this.strm.avail_out = 0;
                var n = s.deflateInit2(this.strm, e.level, e.method, e.windowBits, e.memLevel, e.strategy);
                if (n !== _)
                    throw new Error(c[n]);
                if (e.header && s.deflateSetHeader(this.strm, e.header),
                    e.dictionary) {
                    var i;
                    if (i = "string" == typeof e.dictionary ? l.string2buf(e.dictionary) : "[object ArrayBuffer]" === h.call(e.dictionary) ? new Uint8Array(e.dictionary) : e.dictionary,
                        n = s.deflateSetDictionary(this.strm, i),
                    n !== _)
                        throw new Error(c[n]);
                    this._dict_set = !0
                }
            }
            function i(t, e) {
                var n = new r(e);
                if (n.push(t, !0),
                    n.err)
                    throw n.msg || c[n.err];
                return n.result
            }
            function o(t, e) {
                return e = e || {},
                    e.raw = !0,
                    i(t, e)
            }
            function a(t, e) {
                return e = e || {},
                    e.gzip = !0,
                    i(t, e)
            }
            var s = n(13)
                , u = n(14)
                , l = n(19)
                , c = n(18)
                , f = n(20)
                , h = Object.prototype.toString
                , d = 0
                , p = 4
                , _ = 0
                , y = 1
                , b = 2
                , g = -1
                , v = 0
                , w = 8;
            r.prototype.push = function(t, e) {
                var n, r, i = this.strm, o = this.options.chunkSize;
                if (this.ended)
                    return !1;
                r = e === ~~e ? e : e === !0 ? p : d,
                    "string" == typeof t ? i.input = l.string2buf(t) : "[object ArrayBuffer]" === h.call(t) ? i.input = new Uint8Array(t) : i.input = t,
                    i.next_in = 0,
                    i.avail_in = i.input.length;
                do {
                    if (0 === i.avail_out && (i.output = new u.Buf8(o),
                        i.next_out = 0,
                        i.avail_out = o),
                        n = s.deflate(i, r),
                    n !== y && n !== _)
                        return this.onEnd(n),
                            this.ended = !0,
                            !1;
                    0 !== i.avail_out && (0 !== i.avail_in || r !== p && r !== b) || ("string" === this.options.to ? this.onData(l.buf2binstring(u.shrinkBuf(i.output, i.next_out))) : this.onData(u.shrinkBuf(i.output, i.next_out)))
                } while ((i.avail_in > 0 || 0 === i.avail_out) && n !== y);return r === p ? (n = s.deflateEnd(this.strm),
                    this.onEnd(n),
                    this.ended = !0,
                n === _) : r !== b || (this.onEnd(_),
                    i.avail_out = 0,
                    !0)
            }
                ,
                r.prototype.onData = function(t) {
                    this.chunks.push(t)
                }
                ,
                r.prototype.onEnd = function(t) {
                    t === _ && ("string" === this.options.to ? this.result = this.chunks.join("") : this.result = u.flattenChunks(this.chunks)),
                        this.chunks = [],
                        this.err = t,
                        this.msg = this.strm.msg
                }
                ,
                e.Deflate = r,
                e.deflate = i,
                e.deflateRaw = o,
                e.gzip = a
        }
        , function(t, e, n) {
            "use strict";
            function r(t, e) {
                return t.msg = R[e],
                    e
            }
            function i(t) {
                return (t << 1) - (t > 4 ? 9 : 0)
            }
            function o(t) {
                for (var e = t.length; --e >= 0; )
                    t[e] = 0
            }
            function a(t) {
                var e = t.state
                    , n = e.pending;
                n > t.avail_out && (n = t.avail_out),
                0 !== n && (C.arraySet(t.output, e.pending_buf, e.pending_out, n, t.next_out),
                    t.next_out += n,
                    e.pending_out += n,
                    t.total_out += n,
                    t.avail_out -= n,
                    e.pending -= n,
                0 === e.pending && (e.pending_out = 0))
            }
            function s(t, e) {
                S._tr_flush_block(t, t.block_start >= 0 ? t.block_start : -1, t.strstart - t.block_start, e),
                    t.block_start = t.strstart,
                    a(t.strm)
            }
            function u(t, e) {
                t.pending_buf[t.pending++] = e
            }
            function l(t, e) {
                t.pending_buf[t.pending++] = e >>> 8 & 255,
                    t.pending_buf[t.pending++] = 255 & e
            }
            function c(t, e, n, r) {
                var i = t.avail_in;
                return i > r && (i = r),
                    0 === i ? 0 : (t.avail_in -= i,
                        C.arraySet(e, t.input, t.next_in, i, n),
                        1 === t.state.wrap ? t.adler = z(t.adler, e, i, n) : 2 === t.state.wrap && (t.adler = N(t.adler, e, i, n)),
                        t.next_in += i,
                        t.total_in += i,
                        i)
            }
            function f(t, e) {
                var n, r, i = t.max_chain_length, o = t.strstart, a = t.prev_length, s = t.nice_match, u = t.strstart > t.w_size - ft ? t.strstart - (t.w_size - ft) : 0, l = t.window, c = t.w_mask, f = t.prev, h = t.strstart + ct, d = l[o + a - 1], p = l[o + a];
                t.prev_length >= t.good_match && (i >>= 2),
                s > t.lookahead && (s = t.lookahead);
                do
                    if (n = e,
                    l[n + a] === p && l[n + a - 1] === d && l[n] === l[o] && l[++n] === l[o + 1]) {
                        o += 2,
                            n++;
                        do
                            ;
                        while (l[++o] === l[++n] && l[++o] === l[++n] && l[++o] === l[++n] && l[++o] === l[++n] && l[++o] === l[++n] && l[++o] === l[++n] && l[++o] === l[++n] && l[++o] === l[++n] && o < h);if (r = ct - (h - o),
                            o = h - ct,
                        r > a) {
                            if (t.match_start = e,
                                a = r,
                            r >= s)
                                break;
                            d = l[o + a - 1],
                                p = l[o + a]
                        }
                    }
                while ((e = f[e & c]) > u && 0 !== --i);return a <= t.lookahead ? a : t.lookahead
            }
            function h(t) {
                var e, n, r, i, o, a = t.w_size;
                do {
                    if (i = t.window_size - t.lookahead - t.strstart,
                    t.strstart >= a + (a - ft)) {
                        C.arraySet(t.window, t.window, a, a, 0),
                            t.match_start -= a,
                            t.strstart -= a,
                            t.block_start -= a,
                            n = t.hash_size,
                            e = n;
                        do
                            r = t.head[--e],
                                t.head[e] = r >= a ? r - a : 0;
                        while (--n);n = a,
                            e = n;
                        do
                            r = t.prev[--e],
                                t.prev[e] = r >= a ? r - a : 0;
                        while (--n);i += a
                    }
                    if (0 === t.strm.avail_in)
                        break;
                    if (n = c(t.strm, t.window, t.strstart + t.lookahead, i),
                        t.lookahead += n,
                    t.lookahead + t.insert >= lt)
                        for (o = t.strstart - t.insert,
                                t.ins_h = t.window[o],
                                t.ins_h = (t.ins_h << t.hash_shift ^ t.window[o + 1]) & t.hash_mask; t.insert && (t.ins_h = (t.ins_h << t.hash_shift ^ t.window[o + lt - 1]) & t.hash_mask,
                            t.prev[o & t.w_mask] = t.head[t.ins_h],
                            t.head[t.ins_h] = o,
                            o++,
                            t.insert--,
                            !(t.lookahead + t.insert < lt)); )
                            ;
                } while (t.lookahead < ft && 0 !== t.strm.avail_in)
            }
            function d(t, e) {
                var n = 65535;
                for (n > t.pending_buf_size - 5 && (n = t.pending_buf_size - 5); ; ) {
                    if (t.lookahead <= 1) {
                        if (h(t),
                        0 === t.lookahead && e === I)
                            return wt;
                        if (0 === t.lookahead)
                            break
                    }
                    t.strstart += t.lookahead,
                        t.lookahead = 0;
                    var r = t.block_start + n;
                    if ((0 === t.strstart || t.strstart >= r) && (t.lookahead = t.strstart - r,
                        t.strstart = r,
                        s(t, !1),
                    0 === t.strm.avail_out))
                        return wt;
                    if (t.strstart - t.block_start >= t.w_size - ft && (s(t, !1),
                    0 === t.strm.avail_out))
                        return wt
                }
                return t.insert = 0,
                    e === M ? (s(t, !0),
                        0 === t.strm.avail_out ? Et : kt) : t.strstart > t.block_start && (s(t, !1),
                    0 === t.strm.avail_out) ? wt : wt
            }
            function p(t, e) {
                for (var n, r; ; ) {
                    if (t.lookahead < ft) {
                        if (h(t),
                        t.lookahead < ft && e === I)
                            return wt;
                        if (0 === t.lookahead)
                            break
                    }
                    if (n = 0,
                    t.lookahead >= lt && (t.ins_h = (t.ins_h << t.hash_shift ^ t.window[t.strstart + lt - 1]) & t.hash_mask,
                        n = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h],
                        t.head[t.ins_h] = t.strstart),
                    0 !== n && t.strstart - n <= t.w_size - ft && (t.match_length = f(t, n)),
                    t.match_length >= lt)
                        if (r = S._tr_tally(t, t.strstart - t.match_start, t.match_length - lt),
                            t.lookahead -= t.match_length,
                        t.match_length <= t.max_lazy_match && t.lookahead >= lt) {
                            t.match_length--;
                            do
                                t.strstart++,
                                    t.ins_h = (t.ins_h << t.hash_shift ^ t.window[t.strstart + lt - 1]) & t.hash_mask,
                                    n = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h],
                                    t.head[t.ins_h] = t.strstart;
                            while (0 !== --t.match_length);t.strstart++
                        } else
                            t.strstart += t.match_length,
                                t.match_length = 0,
                                t.ins_h = t.window[t.strstart],
                                t.ins_h = (t.ins_h << t.hash_shift ^ t.window[t.strstart + 1]) & t.hash_mask;
                    else
                        r = S._tr_tally(t, 0, t.window[t.strstart]),
                            t.lookahead--,
                            t.strstart++;
                    if (r && (s(t, !1),
                    0 === t.strm.avail_out))
                        return wt
                }
                return t.insert = t.strstart < lt - 1 ? t.strstart : lt - 1,
                    e === M ? (s(t, !0),
                        0 === t.strm.avail_out ? Et : kt) : t.last_lit && (s(t, !1),
                    0 === t.strm.avail_out) ? wt : mt
            }
            function _(t, e) {
                for (var n, r, i; ; ) {
                    if (t.lookahead < ft) {
                        if (h(t),
                        t.lookahead < ft && e === I)
                            return wt;
                        if (0 === t.lookahead)
                            break
                    }
                    if (n = 0,
                    t.lookahead >= lt && (t.ins_h = (t.ins_h << t.hash_shift ^ t.window[t.strstart + lt - 1]) & t.hash_mask,
                        n = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h],
                        t.head[t.ins_h] = t.strstart),
                        t.prev_length = t.match_length,
                        t.prev_match = t.match_start,
                        t.match_length = lt - 1,
                    0 !== n && t.prev_length < t.max_lazy_match && t.strstart - n <= t.w_size - ft && (t.match_length = f(t, n),
                    t.match_length <= 5 && (t.strategy === V || t.match_length === lt && t.strstart - t.match_start > 4096) && (t.match_length = lt - 1)),
                    t.prev_length >= lt && t.match_length <= t.prev_length) {
                        i = t.strstart + t.lookahead - lt,
                            r = S._tr_tally(t, t.strstart - 1 - t.prev_match, t.prev_length - lt),
                            t.lookahead -= t.prev_length - 1,
                            t.prev_length -= 2;
                        do
                            ++t.strstart <= i && (t.ins_h = (t.ins_h << t.hash_shift ^ t.window[t.strstart + lt - 1]) & t.hash_mask,
                                n = t.prev[t.strstart & t.w_mask] = t.head[t.ins_h],
                                t.head[t.ins_h] = t.strstart);
                        while (0 !== --t.prev_length);if (t.match_available = 0,
                            t.match_length = lt - 1,
                            t.strstart++,
                        r && (s(t, !1),
                        0 === t.strm.avail_out))
                            return wt
                    } else if (t.match_available) {
                        if (r = S._tr_tally(t, 0, t.window[t.strstart - 1]),
                        r && s(t, !1),
                            t.strstart++,
                            t.lookahead--,
                        0 === t.strm.avail_out)
                            return wt
                    } else
                        t.match_available = 1,
                            t.strstart++,
                            t.lookahead--
                }
                return t.match_available && (r = S._tr_tally(t, 0, t.window[t.strstart - 1]),
                    t.match_available = 0),
                    t.insert = t.strstart < lt - 1 ? t.strstart : lt - 1,
                    e === M ? (s(t, !0),
                        0 === t.strm.avail_out ? Et : kt) : t.last_lit && (s(t, !1),
                    0 === t.strm.avail_out) ? wt : mt
            }
            function y(t, e) {
                for (var n, r, i, o, a = t.window; ; ) {
                    if (t.lookahead <= ct) {
                        if (h(t),
                        t.lookahead <= ct && e === I)
                            return wt;
                        if (0 === t.lookahead)
                            break
                    }
                    if (t.match_length = 0,
                    t.lookahead >= lt && t.strstart > 0 && (i = t.strstart - 1,
                        r = a[i],
                    r === a[++i] && r === a[++i] && r === a[++i])) {
                        o = t.strstart + ct;
                        do
                            ;
                        while (r === a[++i] && r === a[++i] && r === a[++i] && r === a[++i] && r === a[++i] && r === a[++i] && r === a[++i] && r === a[++i] && i < o);t.match_length = ct - (o - i),
                        t.match_length > t.lookahead && (t.match_length = t.lookahead)
                    }
                    if (t.match_length >= lt ? (n = S._tr_tally(t, 1, t.match_length - lt),
                        t.lookahead -= t.match_length,
                        t.strstart += t.match_length,
                        t.match_length = 0) : (n = S._tr_tally(t, 0, t.window[t.strstart]),
                        t.lookahead--,
                        t.strstart++),
                    n && (s(t, !1),
                    0 === t.strm.avail_out))
                        return wt
                }
                return t.insert = 0,
                    e === M ? (s(t, !0),
                        0 === t.strm.avail_out ? Et : kt) : t.last_lit && (s(t, !1),
                    0 === t.strm.avail_out) ? wt : mt
            }
            function b(t, e) {
                for (var n; ; ) {
                    if (0 === t.lookahead && (h(t),
                    0 === t.lookahead)) {
                        if (e === I)
                            return wt;
                        break
                    }
                    if (t.match_length = 0,
                        n = S._tr_tally(t, 0, t.window[t.strstart]),
                        t.lookahead--,
                        t.strstart++,
                    n && (s(t, !1),
                    0 === t.strm.avail_out))
                        return wt
                }
                return t.insert = 0,
                    e === M ? (s(t, !0),
                        0 === t.strm.avail_out ? Et : kt) : t.last_lit && (s(t, !1),
                    0 === t.strm.avail_out) ? wt : mt
            }
            function g(t, e, n, r, i) {
                this.good_length = t,
                    this.max_lazy = e,
                    this.nice_length = n,
                    this.max_chain = r,
                    this.func = i
            }
            function v(t) {
                t.window_size = 2 * t.w_size,
                    o(t.head),
                    t.max_lazy_match = P[t.level].max_lazy,
                    t.good_match = P[t.level].good_length,
                    t.nice_match = P[t.level].nice_length,
                    t.max_chain_length = P[t.level].max_chain,
                    t.strstart = 0,
                    t.block_start = 0,
                    t.lookahead = 0,
                    t.insert = 0,
                    t.match_length = t.prev_length = lt - 1,
                    t.match_available = 0,
                    t.ins_h = 0
            }
            function w() {
                this.strm = null,
                    this.status = 0,
                    this.pending_buf = null,
                    this.pending_buf_size = 0,
                    this.pending_out = 0,
                    this.pending = 0,
                    this.wrap = 0,
                    this.gzhead = null,
                    this.gzindex = 0,
                    this.method = q,
                    this.last_flush = -1,
                    this.w_size = 0,
                    this.w_bits = 0,
                    this.w_mask = 0,
                    this.window = null,
                    this.window_size = 0,
                    this.prev = null,
                    this.head = null,
                    this.ins_h = 0,
                    this.hash_size = 0,
                    this.hash_bits = 0,
                    this.hash_mask = 0,
                    this.hash_shift = 0,
                    this.block_start = 0,
                    this.match_length = 0,
                    this.prev_match = 0,
                    this.match_available = 0,
                    this.strstart = 0,
                    this.match_start = 0,
                    this.lookahead = 0,
                    this.prev_length = 0,
                    this.max_chain_length = 0,
                    this.max_lazy_match = 0,
                    this.level = 0,
                    this.strategy = 0,
                    this.good_match = 0,
                    this.nice_match = 0,
                    this.dyn_ltree = new C.Buf16(2 * st),
                    this.dyn_dtree = new C.Buf16(2 * (2 * ot + 1)),
                    this.bl_tree = new C.Buf16(2 * (2 * at + 1)),
                    o(this.dyn_ltree),
                    o(this.dyn_dtree),
                    o(this.bl_tree),
                    this.l_desc = null,
                    this.d_desc = null,
                    this.bl_desc = null,
                    this.bl_count = new C.Buf16(ut + 1),
                    this.heap = new C.Buf16(2 * it + 1),
                    o(this.heap),
                    this.heap_len = 0,
                    this.heap_max = 0,
                    this.depth = new C.Buf16(2 * it + 1),
                    o(this.depth),
                    this.l_buf = 0,
                    this.lit_bufsize = 0,
                    this.last_lit = 0,
                    this.d_buf = 0,
                    this.opt_len = 0,
                    this.static_len = 0,
                    this.matches = 0,
                    this.insert = 0,
                    this.bi_buf = 0,
                    this.bi_valid = 0
            }
            function m(t) {
                var e;
                return t && t.state ? (t.total_in = t.total_out = 0,
                    t.data_type = Q,
                    e = t.state,
                    e.pending = 0,
                    e.pending_out = 0,
                e.wrap < 0 && (e.wrap = -e.wrap),
                    e.status = e.wrap ? dt : gt,
                    t.adler = 2 === e.wrap ? 0 : 1,
                    e.last_flush = I,
                    S._tr_init(e),
                    Y) : r(t, D)
            }
            function E(t) {
                var e = m(t);
                return e === Y && v(t.state),
                    e
            }
            function k(t, e) {
                return t && t.state ? 2 !== t.state.wrap ? D : (t.state.gzhead = e,
                    Y) : D
            }
            function O(t, e, n, i, o, a) {
                if (!t)
                    return D;
                var s = 1;
                if (e === H && (e = 6),
                    i < 0 ? (s = 0,
                        i = -i) : i > 15 && (s = 2,
                        i -= 16),
                o < 1 || o > K || n !== q || i < 8 || i > 15 || e < 0 || e > 9 || a < 0 || a > G)
                    return r(t, D);
                8 === i && (i = 9);
                var u = new w;
                return t.state = u,
                    u.strm = t,
                    u.wrap = s,
                    u.gzhead = null,
                    u.w_bits = i,
                    u.w_size = 1 << u.w_bits,
                    u.w_mask = u.w_size - 1,
                    u.hash_bits = o + 7,
                    u.hash_size = 1 << u.hash_bits,
                    u.hash_mask = u.hash_size - 1,
                    u.hash_shift = ~~((u.hash_bits + lt - 1) / lt),
                    u.window = new C.Buf8(2 * u.w_size),
                    u.head = new C.Buf16(u.hash_size),
                    u.prev = new C.Buf16(u.w_size),
                    u.lit_bufsize = 1 << o + 6,
                    u.pending_buf_size = 4 * u.lit_bufsize,
                    u.pending_buf = new C.Buf8(u.pending_buf_size),
                    u.d_buf = 1 * u.lit_bufsize,
                    u.l_buf = 3 * u.lit_bufsize,
                    u.level = e,
                    u.strategy = a,
                    u.method = n,
                    E(t)
            }
            function T(t, e) {
                return O(t, e, q, tt, et, $)
            }
            function j(t, e) {
                var n, s, c, f;
                if (!t || !t.state || e > F || e < 0)
                    return t ? r(t, D) : D;
                if (s = t.state,
                !t.output || !t.input && 0 !== t.avail_in || s.status === vt && e !== M)
                    return r(t, 0 === t.avail_out ? Z : D);
                if (s.strm = t,
                    n = s.last_flush,
                    s.last_flush = e,
                s.status === dt)
                    if (2 === s.wrap)
                        t.adler = 0,
                            u(s, 31),
                            u(s, 139),
                            u(s, 8),
                            s.gzhead ? (u(s, (s.gzhead.text ? 1 : 0) + (s.gzhead.hcrc ? 2 : 0) + (s.gzhead.extra ? 4 : 0) + (s.gzhead.name ? 8 : 0) + (s.gzhead.comment ? 16 : 0)),
                                u(s, 255 & s.gzhead.time),
                                u(s, s.gzhead.time >> 8 & 255),
                                u(s, s.gzhead.time >> 16 & 255),
                                u(s, s.gzhead.time >> 24 & 255),
                                u(s, 9 === s.level ? 2 : s.strategy >= J || s.level < 2 ? 4 : 0),
                                u(s, 255 & s.gzhead.os),
                            s.gzhead.extra && s.gzhead.extra.length && (u(s, 255 & s.gzhead.extra.length),
                                u(s, s.gzhead.extra.length >> 8 & 255)),
                            s.gzhead.hcrc && (t.adler = N(t.adler, s.pending_buf, s.pending, 0)),
                                s.gzindex = 0,
                                s.status = pt) : (u(s, 0),
                                u(s, 0),
                                u(s, 0),
                                u(s, 0),
                                u(s, 0),
                                u(s, 9 === s.level ? 2 : s.strategy >= J || s.level < 2 ? 4 : 0),
                                u(s, Ot),
                                s.status = gt);
                    else {
                        var h = q + (s.w_bits - 8 << 4) << 8
                            , d = -1;
                        d = s.strategy >= J || s.level < 2 ? 0 : s.level < 6 ? 1 : 6 === s.level ? 2 : 3,
                            h |= d << 6,
                        0 !== s.strstart && (h |= ht),
                            h += 31 - h % 31,
                            s.status = gt,
                            l(s, h),
                        0 !== s.strstart && (l(s, t.adler >>> 16),
                            l(s, 65535 & t.adler)),
                            t.adler = 1
                    }
                if (s.status === pt)
                    if (s.gzhead.extra) {
                        for (c = s.pending; s.gzindex < (65535 & s.gzhead.extra.length) && (s.pending !== s.pending_buf_size || (s.gzhead.hcrc && s.pending > c && (t.adler = N(t.adler, s.pending_buf, s.pending - c, c)),
                            a(t),
                            c = s.pending,
                        s.pending !== s.pending_buf_size)); )
                            u(s, 255 & s.gzhead.extra[s.gzindex]),
                                s.gzindex++;
                        s.gzhead.hcrc && s.pending > c && (t.adler = N(t.adler, s.pending_buf, s.pending - c, c)),
                        s.gzindex === s.gzhead.extra.length && (s.gzindex = 0,
                            s.status = _t)
                    } else
                        s.status = _t;
                if (s.status === _t)
                    if (s.gzhead.name) {
                        c = s.pending;
                        do {
                            if (s.pending === s.pending_buf_size && (s.gzhead.hcrc && s.pending > c && (t.adler = N(t.adler, s.pending_buf, s.pending - c, c)),
                                a(t),
                                c = s.pending,
                            s.pending === s.pending_buf_size)) {
                                f = 1;
                                break
                            }
                            f = s.gzindex < s.gzhead.name.length ? 255 & s.gzhead.name.charCodeAt(s.gzindex++) : 0,
                                u(s, f)
                        } while (0 !== f);s.gzhead.hcrc && s.pending > c && (t.adler = N(t.adler, s.pending_buf, s.pending - c, c)),
                        0 === f && (s.gzindex = 0,
                            s.status = yt)
                    } else
                        s.status = yt;
                if (s.status === yt)
                    if (s.gzhead.comment) {
                        c = s.pending;
                        do {
                            if (s.pending === s.pending_buf_size && (s.gzhead.hcrc && s.pending > c && (t.adler = N(t.adler, s.pending_buf, s.pending - c, c)),
                                a(t),
                                c = s.pending,
                            s.pending === s.pending_buf_size)) {
                                f = 1;
                                break
                            }
                            f = s.gzindex < s.gzhead.comment.length ? 255 & s.gzhead.comment.charCodeAt(s.gzindex++) : 0,
                                u(s, f)
                        } while (0 !== f);s.gzhead.hcrc && s.pending > c && (t.adler = N(t.adler, s.pending_buf, s.pending - c, c)),
                        0 === f && (s.status = bt)
                    } else
                        s.status = bt;
                if (s.status === bt && (s.gzhead.hcrc ? (s.pending + 2 > s.pending_buf_size && a(t),
                s.pending + 2 <= s.pending_buf_size && (u(s, 255 & t.adler),
                    u(s, t.adler >> 8 & 255),
                    t.adler = 0,
                    s.status = gt)) : s.status = gt),
                0 !== s.pending) {
                    if (a(t),
                    0 === t.avail_out)
                        return s.last_flush = -1,
                            Y
                } else if (0 === t.avail_in && i(e) <= i(n) && e !== M)
                    return r(t, Z);
                if (s.status === vt && 0 !== t.avail_in)
                    return r(t, Z);
                if (0 !== t.avail_in || 0 !== s.lookahead || e !== I && s.status !== vt) {
                    var p = s.strategy === J ? b(s, e) : s.strategy === W ? y(s, e) : P[s.level].func(s, e);
                    if (p !== Et && p !== kt || (s.status = vt),
                    p === wt || p === Et)
                        return 0 === t.avail_out && (s.last_flush = -1),
                            Y;
                    if (p === mt && (e === L ? S._tr_align(s) : e !== F && (S._tr_stored_block(s, 0, 0, !1),
                    e === B && (o(s.head),
                    0 === s.lookahead && (s.strstart = 0,
                        s.block_start = 0,
                        s.insert = 0))),
                        a(t),
                    0 === t.avail_out))
                        return s.last_flush = -1,
                            Y
                }
                return e !== M ? Y : s.wrap <= 0 ? U : (2 === s.wrap ? (u(s, 255 & t.adler),
                    u(s, t.adler >> 8 & 255),
                    u(s, t.adler >> 16 & 255),
                    u(s, t.adler >> 24 & 255),
                    u(s, 255 & t.total_in),
                    u(s, t.total_in >> 8 & 255),
                    u(s, t.total_in >> 16 & 255),
                    u(s, t.total_in >> 24 & 255)) : (l(s, t.adler >>> 16),
                    l(s, 65535 & t.adler)),
                    a(t),
                s.wrap > 0 && (s.wrap = -s.wrap),
                    0 !== s.pending ? Y : U)
            }
            function A(t) {
                var e;
                return t && t.state ? (e = t.state.status,
                    e !== dt && e !== pt && e !== _t && e !== yt && e !== bt && e !== gt && e !== vt ? r(t, D) : (t.state = null,
                        e === gt ? r(t, X) : Y)) : D
            }
            function x(t, e) {
                var n, r, i, a, s, u, l, c, f = e.length;
                if (!t || !t.state)
                    return D;
                if (n = t.state,
                    a = n.wrap,
                2 === a || 1 === a && n.status !== dt || n.lookahead)
                    return D;
                for (1 === a && (t.adler = z(t.adler, e, f, 0)),
                        n.wrap = 0,
                    f >= n.w_size && (0 === a && (o(n.head),
                        n.strstart = 0,
                        n.block_start = 0,
                        n.insert = 0),
                        c = new C.Buf8(n.w_size),
                        C.arraySet(c, e, f - n.w_size, n.w_size, 0),
                        e = c,
                        f = n.w_size),
                        s = t.avail_in,
                        u = t.next_in,
                        l = t.input,
                        t.avail_in = f,
                        t.next_in = 0,
                        t.input = e,
                        h(n); n.lookahead >= lt; ) {
                    r = n.strstart,
                        i = n.lookahead - (lt - 1);
                    do
                        n.ins_h = (n.ins_h << n.hash_shift ^ n.window[r + lt - 1]) & n.hash_mask,
                            n.prev[r & n.w_mask] = n.head[n.ins_h],
                            n.head[n.ins_h] = r,
                            r++;
                    while (--i);n.strstart = r,
                        n.lookahead = lt - 1,
                        h(n)
                }
                return n.strstart += n.lookahead,
                    n.block_start = n.strstart,
                    n.insert = n.lookahead,
                    n.lookahead = 0,
                    n.match_length = n.prev_length = lt - 1,
                    n.match_available = 0,
                    t.next_in = u,
                    t.input = l,
                    t.avail_in = s,
                    n.wrap = a,
                    Y
            }
            var P, C = n(14), S = n(15), z = n(16), N = n(17), R = n(18), I = 0, L = 1, B = 3, M = 4, F = 5, Y = 0, U = 1, D = -2, X = -3, Z = -5, H = -1, V = 1, J = 2, W = 3, G = 4, $ = 0, Q = 2, q = 8, K = 9, tt = 15, et = 8, nt = 29, rt = 256, it = rt + 1 + nt, ot = 30, at = 19, st = 2 * it + 1, ut = 15, lt = 3, ct = 258, ft = ct + lt + 1, ht = 32, dt = 42, pt = 69, _t = 73, yt = 91, bt = 103, gt = 113, vt = 666, wt = 1, mt = 2, Et = 3, kt = 4, Ot = 3;
            P = [new g(0,0,0,0,d), new g(4,4,8,4,p), new g(4,5,16,8,p), new g(4,6,32,32,p), new g(4,4,16,16,_), new g(8,16,32,32,_), new g(8,16,128,128,_), new g(8,32,128,256,_), new g(32,128,258,1024,_), new g(32,258,258,4096,_)],
                e.deflateInit = T,
                e.deflateInit2 = O,
                e.deflateReset = E,
                e.deflateResetKeep = m,
                e.deflateSetHeader = k,
                e.deflate = j,
                e.deflateEnd = A,
                e.deflateSetDictionary = x,
                e.deflateInfo = "pako deflate (from Nodeca project)"
        }
        , function(t, e) {
            "use strict";
            var n = "undefined" != typeof Uint8Array && "undefined" != typeof Uint16Array && "undefined" != typeof Int32Array;
            e.assign = function(t) {
                for (var e = Array.prototype.slice.call(arguments, 1); e.length; ) {
                    var n = e.shift();
                    if (n) {
                        if ("object" != typeof n)
                            throw new TypeError(n + "must be non-object");
                        for (var r in n)
                            n.hasOwnProperty(r) && (t[r] = n[r])
                    }
                }
                return t
            }
                ,
                e.shrinkBuf = function(t, e) {
                    return t.length === e ? t : t.subarray ? t.subarray(0, e) : (t.length = e,
                        t)
                }
            ;
            var r = {
                arraySet: function(t, e, n, r, i) {
                    if (e.subarray && t.subarray)
                        return void t.set(e.subarray(n, n + r), i);
                    for (var o = 0; o < r; o++)
                        t[i + o] = e[n + o]
                },
                flattenChunks: function(t) {
                    var e, n, r, i, o, a;
                    for (r = 0,
                            e = 0,
                            n = t.length; e < n; e++)
                        r += t[e].length;
                    for (a = new Uint8Array(r),
                            i = 0,
                            e = 0,
                            n = t.length; e < n; e++)
                        o = t[e],
                            a.set(o, i),
                            i += o.length;
                    return a
                }
            }
                , i = {
                arraySet: function(t, e, n, r, i) {
                    for (var o = 0; o < r; o++)
                        t[i + o] = e[n + o]
                },
                flattenChunks: function(t) {
                    return [].concat.apply([], t)
                }
            };
            e.setTyped = function(t) {
                t ? (e.Buf8 = Uint8Array,
                    e.Buf16 = Uint16Array,
                    e.Buf32 = Int32Array,
                    e.assign(e, r)) : (e.Buf8 = Array,
                    e.Buf16 = Array,
                    e.Buf32 = Array,
                    e.assign(e, i))
            }
                ,
                e.setTyped(n)
        }
        , function(t, e, n) {
            "use strict";
            function r(t) {
                for (var e = t.length; --e >= 0; )
                    t[e] = 0
            }
            function i(t, e, n, r, i) {
                this.static_tree = t,
                    this.extra_bits = e,
                    this.extra_base = n,
                    this.elems = r,
                    this.max_length = i,
                    this.has_stree = t && t.length
            }
            function o(t, e) {
                this.dyn_tree = t,
                    this.max_code = 0,
                    this.stat_desc = e
            }
            function a(t) {
                return t < 256 ? ut[t] : ut[256 + (t >>> 7)]
            }
            function s(t, e) {
                t.pending_buf[t.pending++] = 255 & e,
                    t.pending_buf[t.pending++] = e >>> 8 & 255
            }
            function u(t, e, n) {
                t.bi_valid > G - n ? (t.bi_buf |= e << t.bi_valid & 65535,
                    s(t, t.bi_buf),
                    t.bi_buf = e >> G - t.bi_valid,
                    t.bi_valid += n - G) : (t.bi_buf |= e << t.bi_valid & 65535,
                    t.bi_valid += n)
            }
            function l(t, e, n) {
                u(t, n[2 * e], n[2 * e + 1])
            }
            function c(t, e) {
                var n = 0;
                do
                    n |= 1 & t,
                        t >>>= 1,
                        n <<= 1;
                while (--e > 0);return n >>> 1
            }
            function f(t) {
                16 === t.bi_valid ? (s(t, t.bi_buf),
                    t.bi_buf = 0,
                    t.bi_valid = 0) : t.bi_valid >= 8 && (t.pending_buf[t.pending++] = 255 & t.bi_buf,
                    t.bi_buf >>= 8,
                    t.bi_valid -= 8)
            }
            function h(t, e) {
                var n, r, i, o, a, s, u = e.dyn_tree, l = e.max_code, c = e.stat_desc.static_tree, f = e.stat_desc.has_stree, h = e.stat_desc.extra_bits, d = e.stat_desc.extra_base, p = e.stat_desc.max_length, _ = 0;
                for (o = 0; o <= W; o++)
                    t.bl_count[o] = 0;
                for (u[2 * t.heap[t.heap_max] + 1] = 0,
                        n = t.heap_max + 1; n < J; n++)
                    r = t.heap[n],
                        o = u[2 * u[2 * r + 1] + 1] + 1,
                    o > p && (o = p,
                        _++),
                        u[2 * r + 1] = o,
                    r > l || (t.bl_count[o]++,
                        a = 0,
                    r >= d && (a = h[r - d]),
                        s = u[2 * r],
                        t.opt_len += s * (o + a),
                    f && (t.static_len += s * (c[2 * r + 1] + a)));
                if (0 !== _) {
                    do {
                        for (o = p - 1; 0 === t.bl_count[o]; )
                            o--;
                        t.bl_count[o]--,
                            t.bl_count[o + 1] += 2,
                            t.bl_count[p]--,
                            _ -= 2
                    } while (_ > 0);for (o = p; 0 !== o; o--)
                        for (r = t.bl_count[o]; 0 !== r; )
                            i = t.heap[--n],
                            i > l || (u[2 * i + 1] !== o && (t.opt_len += (o - u[2 * i + 1]) * u[2 * i],
                                u[2 * i + 1] = o),
                                r--)
                }
            }
            function d(t, e, n) {
                var r, i, o = new Array(W + 1), a = 0;
                for (r = 1; r <= W; r++)
                    o[r] = a = a + n[r - 1] << 1;
                for (i = 0; i <= e; i++) {
                    var s = t[2 * i + 1];
                    0 !== s && (t[2 * i] = c(o[s]++, s))
                }
            }
            function p() {
                var t, e, n, r, o, a = new Array(W + 1);
                for (n = 0,
                        r = 0; r < D - 1; r++)
                    for (ct[r] = n,
                            t = 0; t < 1 << et[r]; t++)
                        lt[n++] = r;
                for (lt[n - 1] = r,
                        o = 0,
                        r = 0; r < 16; r++)
                    for (ft[r] = o,
                            t = 0; t < 1 << nt[r]; t++)
                        ut[o++] = r;
                for (o >>= 7; r < H; r++)
                    for (ft[r] = o << 7,
                            t = 0; t < 1 << nt[r] - 7; t++)
                        ut[256 + o++] = r;
                for (e = 0; e <= W; e++)
                    a[e] = 0;
                for (t = 0; t <= 143; )
                    at[2 * t + 1] = 8,
                        t++,
                        a[8]++;
                for (; t <= 255; )
                    at[2 * t + 1] = 9,
                        t++,
                        a[9]++;
                for (; t <= 279; )
                    at[2 * t + 1] = 7,
                        t++,
                        a[7]++;
                for (; t <= 287; )
                    at[2 * t + 1] = 8,
                        t++,
                        a[8]++;
                for (d(at, Z + 1, a),
                        t = 0; t < H; t++)
                    st[2 * t + 1] = 5,
                        st[2 * t] = c(t, 5);
                ht = new i(at,et,X + 1,Z,W),
                    dt = new i(st,nt,0,H,W),
                    pt = new i(new Array(0),rt,0,V,$)
            }
            function _(t) {
                var e;
                for (e = 0; e < Z; e++)
                    t.dyn_ltree[2 * e] = 0;
                for (e = 0; e < H; e++)
                    t.dyn_dtree[2 * e] = 0;
                for (e = 0; e < V; e++)
                    t.bl_tree[2 * e] = 0;
                t.dyn_ltree[2 * Q] = 1,
                    t.opt_len = t.static_len = 0,
                    t.last_lit = t.matches = 0
            }
            function y(t) {
                t.bi_valid > 8 ? s(t, t.bi_buf) : t.bi_valid > 0 && (t.pending_buf[t.pending++] = t.bi_buf),
                    t.bi_buf = 0,
                    t.bi_valid = 0
            }
            function b(t, e, n, r) {
                y(t),
                r && (s(t, n),
                    s(t, ~n)),
                    z.arraySet(t.pending_buf, t.window, e, n, t.pending),
                    t.pending += n
            }
            function g(t, e, n, r) {
                var i = 2 * e
                    , o = 2 * n;
                return t[i] < t[o] || t[i] === t[o] && r[e] <= r[n]
            }
            function v(t, e, n) {
                for (var r = t.heap[n], i = n << 1; i <= t.heap_len && (i < t.heap_len && g(e, t.heap[i + 1], t.heap[i], t.depth) && i++,
                    !g(e, r, t.heap[i], t.depth)); )
                    t.heap[n] = t.heap[i],
                        n = i,
                        i <<= 1;
                t.heap[n] = r
            }
            function w(t, e, n) {
                var r, i, o, s, c = 0;
                if (0 !== t.last_lit)
                    do
                        r = t.pending_buf[t.d_buf + 2 * c] << 8 | t.pending_buf[t.d_buf + 2 * c + 1],
                            i = t.pending_buf[t.l_buf + c],
                            c++,
                            0 === r ? l(t, i, e) : (o = lt[i],
                                l(t, o + X + 1, e),
                                s = et[o],
                            0 !== s && (i -= ct[o],
                                u(t, i, s)),
                                r--,
                                o = a(r),
                                l(t, o, n),
                                s = nt[o],
                            0 !== s && (r -= ft[o],
                                u(t, r, s)));
                    while (c < t.last_lit);l(t, Q, e)
            }
            function m(t, e) {
                var n, r, i, o = e.dyn_tree, a = e.stat_desc.static_tree, s = e.stat_desc.has_stree, u = e.stat_desc.elems, l = -1;
                for (t.heap_len = 0,
                        t.heap_max = J,
                        n = 0; n < u; n++)
                    0 !== o[2 * n] ? (t.heap[++t.heap_len] = l = n,
                        t.depth[n] = 0) : o[2 * n + 1] = 0;
                for (; t.heap_len < 2; )
                    i = t.heap[++t.heap_len] = l < 2 ? ++l : 0,
                        o[2 * i] = 1,
                        t.depth[i] = 0,
                        t.opt_len--,
                    s && (t.static_len -= a[2 * i + 1]);
                for (e.max_code = l,
                        n = t.heap_len >> 1; n >= 1; n--)
                    v(t, o, n);
                i = u;
                do
                    n = t.heap[1],
                        t.heap[1] = t.heap[t.heap_len--],
                        v(t, o, 1),
                        r = t.heap[1],
                        t.heap[--t.heap_max] = n,
                        t.heap[--t.heap_max] = r,
                        o[2 * i] = o[2 * n] + o[2 * r],
                        t.depth[i] = (t.depth[n] >= t.depth[r] ? t.depth[n] : t.depth[r]) + 1,
                        o[2 * n + 1] = o[2 * r + 1] = i,
                        t.heap[1] = i++,
                        v(t, o, 1);
                while (t.heap_len >= 2);t.heap[--t.heap_max] = t.heap[1],
                    h(t, e),
                    d(o, l, t.bl_count)
            }
            function E(t, e, n) {
                var r, i, o = -1, a = e[1], s = 0, u = 7, l = 4;
                for (0 === a && (u = 138,
                    l = 3),
                        e[2 * (n + 1) + 1] = 65535,
                        r = 0; r <= n; r++)
                    i = a,
                        a = e[2 * (r + 1) + 1],
                    ++s < u && i === a || (s < l ? t.bl_tree[2 * i] += s : 0 !== i ? (i !== o && t.bl_tree[2 * i]++,
                        t.bl_tree[2 * q]++) : s <= 10 ? t.bl_tree[2 * K]++ : t.bl_tree[2 * tt]++,
                        s = 0,
                        o = i,
                        0 === a ? (u = 138,
                            l = 3) : i === a ? (u = 6,
                            l = 3) : (u = 7,
                            l = 4))
            }
            function k(t, e, n) {
                var r, i, o = -1, a = e[1], s = 0, c = 7, f = 4;
                for (0 === a && (c = 138,
                    f = 3),
                        r = 0; r <= n; r++)
                    if (i = a,
                        a = e[2 * (r + 1) + 1],
                        !(++s < c && i === a)) {
                        if (s < f) {
                            do
                                l(t, i, t.bl_tree);
                            while (0 !== --s)
                        } else
                            0 !== i ? (i !== o && (l(t, i, t.bl_tree),
                                s--),
                                l(t, q, t.bl_tree),
                                u(t, s - 3, 2)) : s <= 10 ? (l(t, K, t.bl_tree),
                                u(t, s - 3, 3)) : (l(t, tt, t.bl_tree),
                                u(t, s - 11, 7));
                        s = 0,
                            o = i,
                            0 === a ? (c = 138,
                                f = 3) : i === a ? (c = 6,
                                f = 3) : (c = 7,
                                f = 4)
                    }
            }
            function O(t) {
                var e;
                for (E(t, t.dyn_ltree, t.l_desc.max_code),
                        E(t, t.dyn_dtree, t.d_desc.max_code),
                        m(t, t.bl_desc),
                        e = V - 1; e >= 3 && 0 === t.bl_tree[2 * it[e] + 1]; e--)
                    ;
                return t.opt_len += 3 * (e + 1) + 5 + 5 + 4,
                    e
            }
            function T(t, e, n, r) {
                var i;
                for (u(t, e - 257, 5),
                        u(t, n - 1, 5),
                        u(t, r - 4, 4),
                        i = 0; i < r; i++)
                    u(t, t.bl_tree[2 * it[i] + 1], 3);
                k(t, t.dyn_ltree, e - 1),
                    k(t, t.dyn_dtree, n - 1)
            }
            function j(t) {
                var e, n = 4093624447;
                for (e = 0; e <= 31; e++,
                    n >>>= 1)
                    if (1 & n && 0 !== t.dyn_ltree[2 * e])
                        return R;
                if (0 !== t.dyn_ltree[18] || 0 !== t.dyn_ltree[20] || 0 !== t.dyn_ltree[26])
                    return I;
                for (e = 32; e < X; e++)
                    if (0 !== t.dyn_ltree[2 * e])
                        return I;
                return R
            }
            function A(t) {
                _t || (p(),
                    _t = !0),
                    t.l_desc = new o(t.dyn_ltree,ht),
                    t.d_desc = new o(t.dyn_dtree,dt),
                    t.bl_desc = new o(t.bl_tree,pt),
                    t.bi_buf = 0,
                    t.bi_valid = 0,
                    _(t)
            }
            function x(t, e, n, r) {
                u(t, (B << 1) + (r ? 1 : 0), 3),
                    b(t, e, n, !0)
            }
            function P(t) {
                u(t, M << 1, 3),
                    l(t, Q, at),
                    f(t)
            }
            function C(t, e, n, r) {
                var i, o, a = 0;
                t.level > 0 ? (t.strm.data_type === L && (t.strm.data_type = j(t)),
                    m(t, t.l_desc),
                    m(t, t.d_desc),
                    a = O(t),
                    i = t.opt_len + 3 + 7 >>> 3,
                    o = t.static_len + 3 + 7 >>> 3,
                o <= i && (i = o)) : i = o = n + 5,
                    n + 4 <= i && e !== -1 ? x(t, e, n, r) : t.strategy === N || o === i ? (u(t, (M << 1) + (r ? 1 : 0), 3),
                        w(t, at, st)) : (u(t, (F << 1) + (r ? 1 : 0), 3),
                        T(t, t.l_desc.max_code + 1, t.d_desc.max_code + 1, a + 1),
                        w(t, t.dyn_ltree, t.dyn_dtree)),
                    _(t),
                r && y(t)
            }
            function S(t, e, n) {
                return t.pending_buf[t.d_buf + 2 * t.last_lit] = e >>> 8 & 255,
                    t.pending_buf[t.d_buf + 2 * t.last_lit + 1] = 255 & e,
                    t.pending_buf[t.l_buf + t.last_lit] = 255 & n,
                    t.last_lit++,
                    0 === e ? t.dyn_ltree[2 * n]++ : (t.matches++,
                        e--,
                        t.dyn_ltree[2 * (lt[n] + X + 1)]++,
                        t.dyn_dtree[2 * a(e)]++),
                t.last_lit === t.lit_bufsize - 1
            }
            var z = n(14)
                , N = 4
                , R = 0
                , I = 1
                , L = 2
                , B = 0
                , M = 1
                , F = 2
                , Y = 3
                , U = 258
                , D = 29
                , X = 256
                , Z = X + 1 + D
                , H = 30
                , V = 19
                , J = 2 * Z + 1
                , W = 15
                , G = 16
                , $ = 7
                , Q = 256
                , q = 16
                , K = 17
                , tt = 18
                , et = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0]
                , nt = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13]
                , rt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7]
                , it = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
                , ot = 512
                , at = new Array(2 * (Z + 2));
            r(at);
            var st = new Array(2 * H);
            r(st);
            var ut = new Array(ot);
            r(ut);
            var lt = new Array(U - Y + 1);
            r(lt);
            var ct = new Array(D);
            r(ct);
            var ft = new Array(H);
            r(ft);
            var ht, dt, pt, _t = !1;
            e._tr_init = A,
                e._tr_stored_block = x,
                e._tr_flush_block = C,
                e._tr_tally = S,
                e._tr_align = P
        }
        , function(t, e) {
            "use strict";
            function n(t, e, n, r) {
                for (var i = 65535 & t | 0, o = t >>> 16 & 65535 | 0, a = 0; 0 !== n; ) {
                    a = n > 2e3 ? 2e3 : n,
                        n -= a;
                    do
                        i = i + e[r++] | 0,
                            o = o + i | 0;
                    while (--a);i %= 65521,
                        o %= 65521
                }
                return i | o << 16 | 0
            }
            t.exports = n
        }
        , function(t, e) {
            "use strict";
            function n() {
                for (var t, e = [], n = 0; n < 256; n++) {
                    t = n;
                    for (var r = 0; r < 8; r++)
                        t = 1 & t ? 3988292384 ^ t >>> 1 : t >>> 1;
                    e[n] = t
                }
                return e
            }
            function r(t, e, n, r) {
                var o = i
                    , a = r + n;
                t ^= -1;
                for (var s = r; s < a; s++)
                    t = t >>> 8 ^ o[255 & (t ^ e[s])];
                return t ^ -1
            }
            var i = n();
            t.exports = r
        }
        , function(t, e) {
            "use strict";
            t.exports = {
                2: "need dictionary",
                1: "stream end",
                0: "",
                "-1": "file error",
                "-2": "stream error",
                "-3": "data error",
                "-4": "insufficient memory",
                "-5": "buffer error",
                "-6": "incompatible version"
            }
        }
        , function(t, e, n) {
            "use strict";
            function r(t, e) {
                if (e < 65537 && (t.subarray && a || !t.subarray && o))
                    return String.fromCharCode.apply(null, i.shrinkBuf(t, e));
                for (var n = "", r = 0; r < e; r++)
                    n += String.fromCharCode(t[r]);
                return n
            }
            var i = n(14)
                , o = !0
                , a = !0;
            try {
                String.fromCharCode.apply(null, [0])
            } catch (t) {
                o = !1
            }
            try {
                String.fromCharCode.apply(null, new Uint8Array(1))
            } catch (t) {
                a = !1
            }
            for (var s = new i.Buf8(256), u = 0; u < 256; u++)
                s[u] = u >= 252 ? 6 : u >= 248 ? 5 : u >= 240 ? 4 : u >= 224 ? 3 : u >= 192 ? 2 : 1;
            s[254] = s[254] = 1,
                e.string2buf = function(t) {
                    var e, n, r, o, a, s = t.length, u = 0;
                    for (o = 0; o < s; o++)
                        n = t.charCodeAt(o),
                        55296 === (64512 & n) && o + 1 < s && (r = t.charCodeAt(o + 1),
                        56320 === (64512 & r) && (n = 65536 + (n - 55296 << 10) + (r - 56320),
                            o++)),
                            u += n < 128 ? 1 : n < 2048 ? 2 : n < 65536 ? 3 : 4;
                    for (e = new i.Buf8(u),
                            a = 0,
                            o = 0; a < u; o++)
                        n = t.charCodeAt(o),
                        55296 === (64512 & n) && o + 1 < s && (r = t.charCodeAt(o + 1),
                        56320 === (64512 & r) && (n = 65536 + (n - 55296 << 10) + (r - 56320),
                            o++)),
                            n < 128 ? e[a++] = n : n < 2048 ? (e[a++] = 192 | n >>> 6,
                                e[a++] = 128 | 63 & n) : n < 65536 ? (e[a++] = 224 | n >>> 12,
                                e[a++] = 128 | n >>> 6 & 63,
                                e[a++] = 128 | 63 & n) : (e[a++] = 240 | n >>> 18,
                                e[a++] = 128 | n >>> 12 & 63,
                                e[a++] = 128 | n >>> 6 & 63,
                                e[a++] = 128 | 63 & n);
                    return e
                }
                ,
                e.buf2binstring = function(t) {
                    return r(t, t.length)
                }
                ,
                e.binstring2buf = function(t) {
                    for (var e = new i.Buf8(t.length), n = 0, r = e.length; n < r; n++)
                        e[n] = t.charCodeAt(n);
                    return e
                }
                ,
                e.buf2string = function(t, e) {
                    var n, i, o, a, u = e || t.length, l = new Array(2 * u);
                    for (i = 0,
                            n = 0; n < u; )
                        if (o = t[n++],
                        o < 128)
                            l[i++] = o;
                        else if (a = s[o],
                        a > 4)
                            l[i++] = 65533,
                                n += a - 1;
                        else {
                            for (o &= 2 === a ? 31 : 3 === a ? 15 : 7; a > 1 && n < u; )
                                o = o << 6 | 63 & t[n++],
                                    a--;
                            a > 1 ? l[i++] = 65533 : o < 65536 ? l[i++] = o : (o -= 65536,
                                l[i++] = 55296 | o >> 10 & 1023,
                                l[i++] = 56320 | 1023 & o)
                        }
                    return r(l, i)
                }
                ,
                e.utf8border = function(t, e) {
                    var n;
                    for (e = e || t.length,
                        e > t.length && (e = t.length),
                            n = e - 1; n >= 0 && 128 === (192 & t[n]); )
                        n--;
                    return n < 0 ? e : 0 === n ? e : n + s[t[n]] > e ? n : e
                }
        }
        , function(t, e) {
            "use strict";
            function n() {
                this.input = null,
                    this.next_in = 0,
                    this.avail_in = 0,
                    this.total_in = 0,
                    this.output = null,
                    this.next_out = 0,
                    this.avail_out = 0,
                    this.total_out = 0,
                    this.msg = "",
                    this.state = null,
                    this.data_type = 2,
                    this.adler = 0
            }
            t.exports = n
        }
        , function(t, e, n) {
            "use strict";
            e.decode = e.parse = n(22),
                e.encode = e.stringify = n(23)
        }
        , function(t, e) {
            "use strict";
            function n(t, e) {
                return Object.prototype.hasOwnProperty.call(t, e)
            }
            t.exports = function(t, e, r, i) {
                e = e || "&",
                    r = r || "=";
                var o = {};
                if ("string" != typeof t || 0 === t.length)
                    return o;
                var a = /\+/g;
                t = t.split(e);
                var s = 1e3;
                i && "number" == typeof i.maxKeys && (s = i.maxKeys);
                var u = t.length;
                s > 0 && u > s && (u = s);
                for (var l = 0; l < u; ++l) {
                    var c, f, h, d, p = t[l].replace(a, "%20"), _ = p.indexOf(r);
                    _ >= 0 ? (c = p.substr(0, _),
                        f = p.substr(_ + 1)) : (c = p,
                        f = ""),
                        h = decodeURIComponent(c),
                        d = decodeURIComponent(f),
                        n(o, h) ? Array.isArray(o[h]) ? o[h].push(d) : o[h] = [o[h], d] : o[h] = d
                }
                return o
            }
        }
        , function(t, e) {
            "use strict";
            var n = function(t) {
                switch (typeof t) {
                    case "string":
                        return t;
                    case "boolean":
                        return t ? "true" : "false";
                    case "number":
                        return isFinite(t) ? t : "";
                    default:
                        return ""
                }
            };
            t.exports = function(t, e, r, i) {
                return e = e || "&",
                    r = r || "=",
                null === t && (t = void 0),
                    "object" == typeof t ? Object.keys(t).map(function(i) {
                        var o = encodeURIComponent(n(i)) + r;
                        return Array.isArray(t[i]) ? t[i].map(function(t) {
                            return o + encodeURIComponent(n(t))
                        }).join(e) : o + encodeURIComponent(n(t[i]))
                    }).join(e) : i ? encodeURIComponent(n(i)) + r + encodeURIComponent(n(t)) : ""
            }
        }
        , function(t, e) {
            "use strict";
            function n(t, e) {
                return "hasAttribute"in t ? t.hasAttribute(e) : k.filter(t.attributes, function(t) {
                    return t.nodeName === e
                }).length > 0
            }
            function r(t) {
                var e = "d2ViZHJpdmVyLF9fZHJpdmVyX2V2YWx1YXRlLF9fd2ViZHJpdmVyX2V2YWx1YXRlLF9fc2VsZW5pdW1fZXZhbHVhdGUsX19meGRyaXZlcl9ldmFsdWF0ZSxfX2RyaXZlcl91bndyYXBwZWQsX193ZWJkcml2ZXJfdW53cmFwcGVkLF9fc2VsZW5pdW1fdW53cmFwcGVkLF9fZnhkcml2ZXJfdW53cmFwcGVk"
                    , n = k.parse(e);
                return k.filter(n, i(t)).length > 0
            }
            function i(t) {
                return function(e) {
                    return e in t
                }
            }
            function o(t) {
                return atob("X193ZWJkcml2ZXJGdW5j")in t
            }
            function a(t) {
                var e = "d2ViZHJpdmVyLF9TZWxlbml1bV9JREVfUmVjb3JkZXIsX3NlbGVuaXVtLGNhbGxlZFNlbGVuaXVt"
                    , n = k.parse(e);
                return k.filter(n, i(t)).length > 0
            }
            function s(t) {
                return atob("ZG9tQXV0b21hdGlvbg==")in t || atob("ZG9tQXV0b21hdGlvbkNvbnRyb2xsZXI=")in t
            }
            function u(t) {
                return t.documentElement && n(t.documentElement, atob("d2ViZHJpdmVy"))
            }
            function l(t) {
                return atob("X19sYXN0V2F0aXJBbGVydA==")in t || atob("X19sYXN0V2F0aXJDb25maXJt")in t || atob("X19sYXN0V2F0aXJQcm9tcHQ=")in t
            }
            function c(t) {
                return t[atob("d2ViZHJpdmVy")] || !1
            }
            function f(t) {
                return atob("d2ViZHJpdmVy")in t
            }
            function h(t) {
                return atob("X193ZWJkcml2ZXJfc2NyaXB0X2Zu")in t
            }
            function d(t) {
                var e = !1;
                try {
                    e = t.cookie.indexOf(atob("Q2hyb21lRHJpdmVyd2plcnM5MDhmbGpzZGYzNzQ1OWZzZGZnZGZ3cnU9")) > -1
                } catch (t) {}
                return e
            }
            function p(t) {
                return atob("JGNkY19hc2RqZmxhc3V0b3BmaHZjWkxtY2ZsXw==")in t || atob("JGNocm9tZV9hc3luY1NjcmlwdEluZm8=")in t
            }
            function _(t) {
                return atob("X1dFQkRSSVZFUl9FTEVNX0NBQ0hF")in t
            }
            function y(t) {
                return atob("X18kd2ViZHJpdmVyQXN5bmNFeGVjdXRvcg==")in t
            }
            function b(t) {
                var e, n = [];
                for (e = 0; e < t.length; e++)
                    n.push(t[e]);
                return n
            }
            function g(t) {
                return n(t, atob("Y2RfZnJhbWVfaWRf"))
            }
            function v(t) {
                var e = b(t.getElementsByTagName("iframe"))
                    , n = b(t.getElementsByTagName("frame"))
                    , r = e.concat(n)
                    , i = k.filter(r, g);
                return i.length > 0
            }
            function w(t) {
                var e = "ZHJpdmVyLWV2YWx1YXRlLHdlYmRyaXZlci1ldmFsdWF0ZSxzZWxlbml1bS1ldmFsdWF0ZSx3ZWJkcml2ZXJDb21tYW5kLHdlYmRyaXZlci1ldmFsdWF0ZS1yZXNwb25zZQ=="
                    , n = k.parse(e);
                document.addEventListener && k.forEach(n, function(e) {
                    document.addEventListener(e, m(e, t), !1)
                })
            }
            function m(t, e) {
                return function n() {
                    e("lwe"),
                        document.removeEventListener(t, n)
                }
            }
            function E(t) {
                var e = 0
                    , n = setInterval(function() {
                    var r = {};
                    r.f = h(window),
                        r.v = d(document),
                        r.p = p(document),
                        r.h = _(window),
                        r.l = y(document),
                        r.S = v(document);
                    for (var i = k.ownKeys(r), o = 0; o < i.length; o++)
                        if (r[i[o]] === !0) {
                            clearInterval(n),
                                t("lwc" + i[o]);
                            break
                        }
                    ++e > 60 && clearInterval(n)
                }, 500)
            }
            var k = {
                filter: function(t, e) {
                    var n, r = [];
                    for (n = 0; n < t.length; n++)
                        e(t[n], n, t) && r.push(t[n]);
                    return r
                },
                forEach: function(t, e) {
                    var n;
                    for (n = 0; n < t.length; n++)
                        e(t[n], n, t)
                },
                ownKeys: function(t) {
                    var e, n = [];
                    for (e in t)
                        t.hasOwnProperty(e) && n.push(e);
                    return n
                },
                parse: function(t) {
                    return t ? atob(t).split(",") : ""
                }
            }
                , O = function() {
                return u(document) ? "dw" : r(document) ? "de" : a(document) ? "di" : o(window) ? "wf" : s(window) ? "" : l(window) ? "wwt" : f(window) ? "ww" : c(navigator) ? "gw" : ""
            }
                , T = function(t) {
                w(t),
                    E(t)
            };
            t.exports = {
                getwd: O,
                listenwd: T
            }
        }
        , function(t, e, n) {
            "use strict";
            var r = Object.prototype.hasOwnProperty
                , i = Object.prototype.toString
                , o = Array.prototype.slice
                , a = n(26)
                , s = Object.prototype.propertyIsEnumerable
                , u = !s.call({
                toString: null
            }, "toString")
                , l = s.call(function() {}, "prototype")
                , c = ["toString", "toLocaleString", "valueOf", "hasOwnProperty", "isPrototypeOf", "propertyIsEnumerable", "constructor"]
                , f = function(t) {
                var e = t.constructor;
                return e && e.prototype === t
            }
                , h = {
                $console: !0,
                $external: !0,
                $frame: !0,
                $frameElement: !0,
                $frames: !0,
                $innerHeight: !0,
                $innerWidth: !0,
                $outerHeight: !0,
                $outerWidth: !0,
                $pageXOffset: !0,
                $pageYOffset: !0,
                $parent: !0,
                $scrollLeft: !0,
                $scrollTop: !0,
                $scrollX: !0,
                $scrollY: !0,
                $self: !0,
                $webkitIndexedDB: !0,
                $webkitStorageInfo: !0,
                $window: !0
            }
                , d = function() {
                if ("undefined" == typeof window)
                    return !1;
                for (var t in window)
                    try {
                        if (!h["$" + t] && r.call(window, t) && null !== window[t] && "object" == typeof window[t])
                            try {
                                f(window[t])
                            } catch (t) {
                                return !0
                            }
                    } catch (t) {
                        return !0
                    }
                return !1
            }()
                , p = function(t) {
                if ("undefined" == typeof window || !d)
                    return f(t);
                try {
                    return f(t)
                } catch (t) {
                    return !1
                }
            }
                , _ = function(t) {
                var e = null !== t && "object" == typeof t
                    , n = "[object Function]" === i.call(t)
                    , o = a(t)
                    , s = e && "[object String]" === i.call(t)
                    , f = [];
                if (!e && !n && !o)
                    throw new TypeError("Object.keys called on a non-object");
                var h = l && n;
                if (s && t.length > 0 && !r.call(t, 0))
                    for (var d = 0; d < t.length; ++d)
                        f.push(String(d));
                if (o && t.length > 0)
                    for (var _ = 0; _ < t.length; ++_)
                        f.push(String(_));
                else
                    for (var y in t)
                        h && "prototype" === y || !r.call(t, y) || f.push(String(y));
                if (u)
                    for (var b = p(t), g = 0; g < c.length; ++g)
                        b && "constructor" === c[g] || !r.call(t, c[g]) || f.push(c[g]);
                return f
            };
            _.shim = function() {
                if (Object.keys) {
                    var t = function() {
                        return 2 === (Object.keys(arguments) || "").length
                    }(1, 2);
                    if (!t) {
                        var e = Object.keys;
                        Object.keys = function(t) {
                            return e(a(t) ? o.call(t) : t)
                        }
                    }
                } else
                    Object.keys = _;
                return Object.keys || _
            }
                ,
                t.exports = _
        }
        , function(t, e) {
            "use strict";
            var n = Object.prototype.toString;
            t.exports = function(t) {
                var e = n.call(t)
                    , r = "[object Arguments]" === e;
                return r || (r = "[object Array]" !== e && null !== t && "object" == typeof t && "number" == typeof t.length && t.length >= 0 && "[object Function]" === n.call(t.callee)),
                    r
            }
        }
        , function(t, e, n) {
            var r;
            (function(t, i) {
                    (function() {
                            function o(t, e) {
                                function n(t) {
                                    if (n[t] !== y)
                                        return n[t];
                                    var o;
                                    if ("bug-string-char-index" == t)
                                        o = "a" != "a"[0];
                                    else if ("json" == t)
                                        o = n("json-stringify") && n("json-parse");
                                    else {
                                        var a, s = '{"a":[1,true,false,null,"\\u0000\\b\\n\\f\\r\\t"]}';
                                        if ("json-stringify" == t) {
                                            var l = e.stringify
                                                , c = "function" == typeof l && v;
                                            if (c) {
                                                (a = function() {
                                                        return 1
                                                    }
                                                ).toJSON = a;
                                                try {
                                                    c = "0" === l(0) && "0" === l(new r) && '""' == l(new i) && l(g) === y && l(y) === y && l() === y && "1" === l(a) && "[1]" == l([a]) && "[null]" == l([y]) && "null" == l(null) && "[null,null,null]" == l([y, g, null]) && l({
                                                        a: [a, !0, !1, null, "\0\b\n\f\r\t"]
                                                    }) == s && "1" === l(null, a) && "[\n 1,\n 2\n]" == l([1, 2], null, 1) && '"-271821-04-20T00:00:00.000Z"' == l(new u((-864e13))) && '"+275760-09-13T00:00:00.000Z"' == l(new u(864e13)) && '"-000001-01-01T00:00:00.000Z"' == l(new u((-621987552e5))) && '"1969-12-31T23:59:59.999Z"' == l(new u((-1)))
                                                } catch (t) {
                                                    c = !1
                                                }
                                            }
                                            o = c
                                        }
                                        if ("json-parse" == t) {
                                            var f = e.parse;
                                            if ("function" == typeof f)
                                                try {
                                                    if (0 === f("0") && !f(!1)) {
                                                        a = f(s);
                                                        var h = 5 == a.a.length && 1 === a.a[0];
                                                        if (h) {
                                                            try {
                                                                h = !f('"\t"')
                                                            } catch (t) {}
                                                            if (h)
                                                                try {
                                                                    h = 1 !== f("01")
                                                                } catch (t) {}
                                                            if (h)
                                                                try {
                                                                    h = 1 !== f("1.")
                                                                } catch (t) {}
                                                        }
                                                    }
                                                } catch (t) {
                                                    h = !1
                                                }
                                            o = h
                                        }
                                    }
                                    return n[t] = !!o
                                }
                                t || (t = l.Object()),
                                e || (e = l.Object());
                                var r = t.Number || l.Number
                                    , i = t.String || l.String
                                    , a = t.Object || l.Object
                                    , u = t.Date || l.Date
                                    , c = t.SyntaxError || l.SyntaxError
                                    , f = t.TypeError || l.TypeError
                                    , h = t.Math || l.Math
                                    , d = t.JSON || l.JSON;
                                "object" == typeof d && d && (e.stringify = d.stringify,
                                    e.parse = d.parse);
                                var p, _, y, b = a.prototype, g = b.toString, v = new u((-0xc782b5b800cec));
                                try {
                                    v = v.getUTCFullYear() == -109252 && 0 === v.getUTCMonth() && 1 === v.getUTCDate() && 10 == v.getUTCHours() && 37 == v.getUTCMinutes() && 6 == v.getUTCSeconds() && 708 == v.getUTCMilliseconds()
                                } catch (t) {}
                                if (!n("json")) {
                                    var w = "[object Function]"
                                        , m = "[object Date]"
                                        , E = "[object Number]"
                                        , k = "[object String]"
                                        , O = "[object Array]"
                                        , T = "[object Boolean]"
                                        , j = n("bug-string-char-index");
                                    if (!v)
                                        var A = h.floor
                                            , x = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
                                            , P = function(t, e) {
                                            return x[e] + 365 * (t - 1970) + A((t - 1969 + (e = +(e > 1))) / 4) - A((t - 1901 + e) / 100) + A((t - 1601 + e) / 400)
                                        };
                                    if ((p = b.hasOwnProperty) || (p = function(t) {
                                            var e, n = {};
                                            return (n.__proto__ = null,
                                                n.__proto__ = {
                                                    toString: 1
                                                },
                                                n).toString != g ? p = function(t) {
                                                    var e = this.__proto__
                                                        , n = t in (this.__proto__ = null,
                                                        this);
                                                    return this.__proto__ = e,
                                                        n
                                                }
                                                : (e = n.constructor,
                                                        p = function(t) {
                                                            var n = (this.constructor || e).prototype;
                                                            return t in this && !(t in n && this[t] === n[t])
                                                        }
                                                ),
                                                n = null,
                                                p.call(this, t)
                                        }
                                    ),
                                        _ = function(t, e) {
                                            var n, r, i, o = 0;
                                            (n = function() {
                                                    this.valueOf = 0
                                                }
                                            ).prototype.valueOf = 0,
                                                r = new n;
                                            for (i in r)
                                                p.call(r, i) && o++;
                                            return n = r = null,
                                                o ? _ = 2 == o ? function(t, e) {
                                                        var n, r = {}, i = g.call(t) == w;
                                                        for (n in t)
                                                            i && "prototype" == n || p.call(r, n) || !(r[n] = 1) || !p.call(t, n) || e(n)
                                                    }
                                                    : function(t, e) {
                                                        var n, r, i = g.call(t) == w;
                                                        for (n in t)
                                                            i && "prototype" == n || !p.call(t, n) || (r = "constructor" === n) || e(n);
                                                        (r || p.call(t, n = "constructor")) && e(n)
                                                    }
                                                    : (r = ["valueOf", "toString", "toLocaleString", "propertyIsEnumerable", "isPrototypeOf", "hasOwnProperty", "constructor"],
                                                            _ = function(t, e) {
                                                                var n, i, o = g.call(t) == w, a = !o && "function" != typeof t.constructor && s[typeof t.hasOwnProperty] && t.hasOwnProperty || p;
                                                                for (n in t)
                                                                    o && "prototype" == n || !a.call(t, n) || e(n);
                                                                for (i = r.length; n = r[--i]; a.call(t, n) && e(n))
                                                                    ;
                                                            }
                                                    ),
                                                _(t, e)
                                        }
                                        ,
                                        !n("json-stringify")) {
                                        var C = {
                                            92: "\\\\",
                                            34: '\\"',
                                            8: "\\b",
                                            12: "\\f",
                                            10: "\\n",
                                            13: "\\r",
                                            9: "\\t"
                                        }
                                            , S = "000000"
                                            , z = function(t, e) {
                                            return (S + (e || 0)).slice(-t)
                                        }
                                            , N = "\\u00"
                                            , R = function(t) {
                                            for (var e = '"', n = 0, r = t.length, i = !j || r > 10, o = i && (j ? t.split("") : t); n < r; n++) {
                                                var a = t.charCodeAt(n);
                                                switch (a) {
                                                    case 8:
                                                    case 9:
                                                    case 10:
                                                    case 12:
                                                    case 13:
                                                    case 34:
                                                    case 92:
                                                        e += C[a];
                                                        break;
                                                    default:
                                                        if (a < 32) {
                                                            e += N + z(2, a.toString(16));
                                                            break
                                                        }
                                                        e += i ? o[n] : t.charAt(n)
                                                }
                                            }
                                            return e + '"'
                                        }
                                            , I = function(t, e, n, r, i, o, a) {
                                            var s, u, l, c, h, d, b, v, w, j, x, C, S, N, L, B;
                                            try {
                                                s = e[t]
                                            } catch (t) {}
                                            if ("object" == typeof s && s)
                                                if (u = g.call(s),
                                                u != m || p.call(s, "toJSON"))
                                                    "function" == typeof s.toJSON && (u != E && u != k && u != O || p.call(s, "toJSON")) && (s = s.toJSON(t));
                                                else if (s > -1 / 0 && s < 1 / 0) {
                                                    if (P) {
                                                        for (h = A(s / 864e5),
                                                                l = A(h / 365.2425) + 1970 - 1; P(l + 1, 0) <= h; l++)
                                                            ;
                                                        for (c = A((h - P(l, 0)) / 30.42); P(l, c + 1) <= h; c++)
                                                            ;
                                                        h = 1 + h - P(l, c),
                                                            d = (s % 864e5 + 864e5) % 864e5,
                                                            b = A(d / 36e5) % 24,
                                                            v = A(d / 6e4) % 60,
                                                            w = A(d / 1e3) % 60,
                                                            j = d % 1e3
                                                    } else
                                                        l = s.getUTCFullYear(),
                                                            c = s.getUTCMonth(),
                                                            h = s.getUTCDate(),
                                                            b = s.getUTCHours(),
                                                            v = s.getUTCMinutes(),
                                                            w = s.getUTCSeconds(),
                                                            j = s.getUTCMilliseconds();
                                                    s = (l <= 0 || l >= 1e4 ? (l < 0 ? "-" : "+") + z(6, l < 0 ? -l : l) : z(4, l)) + "-" + z(2, c + 1) + "-" + z(2, h) + "T" + z(2, b) + ":" + z(2, v) + ":" + z(2, w) + "." + z(3, j) + "Z"
                                                } else
                                                    s = null;
                                            if (n && (s = n.call(e, t, s)),
                                            null === s)
                                                return "null";
                                            if (u = g.call(s),
                                            u == T)
                                                return "" + s;
                                            if (u == E)
                                                return s > -1 / 0 && s < 1 / 0 ? "" + s : "null";
                                            if (u == k)
                                                return R("" + s);
                                            if ("object" == typeof s) {
                                                for (N = a.length; N--; )
                                                    if (a[N] === s)
                                                        throw f();
                                                if (a.push(s),
                                                    x = [],
                                                    L = o,
                                                    o += i,
                                                u == O) {
                                                    for (S = 0,
                                                            N = s.length; S < N; S++)
                                                        C = I(S, s, n, r, i, o, a),
                                                            x.push(C === y ? "null" : C);
                                                    B = x.length ? i ? "[\n" + o + x.join(",\n" + o) + "\n" + L + "]" : "[" + x.join(",") + "]" : "[]"
                                                } else
                                                    _(r || s, function(t) {
                                                        var e = I(t, s, n, r, i, o, a);
                                                        e !== y && x.push(R(t) + ":" + (i ? " " : "") + e)
                                                    }),
                                                        B = x.length ? i ? "{\n" + o + x.join(",\n" + o) + "\n" + L + "}" : "{" + x.join(",") + "}" : "{}";
                                                return a.pop(),
                                                    B
                                            }
                                        };
                                        e.stringify = function(t, e, n) {
                                            var r, i, o, a;
                                            if (s[typeof e] && e)
                                                if ((a = g.call(e)) == w)
                                                    i = e;
                                                else if (a == O) {
                                                    o = {};
                                                    for (var u, l = 0, c = e.length; l < c; u = e[l++],
                                                        a = g.call(u),
                                                    (a == k || a == E) && (o[u] = 1))
                                                        ;
                                                }
                                            if (n)
                                                if ((a = g.call(n)) == E) {
                                                    if ((n -= n % 1) > 0)
                                                        for (r = "",
                                                            n > 10 && (n = 10); r.length < n; r += " ")
                                                            ;
                                                } else
                                                    a == k && (r = n.length <= 10 ? n : n.slice(0, 10));
                                            return I("", (u = {},
                                                u[""] = t,
                                                u), i, o, r, "", [])
                                        }
                                    }
                                    if (!n("json-parse")) {
                                        var L, B, M = i.fromCharCode, F = {
                                            92: "\\",
                                            34: '"',
                                            47: "/",
                                            98: "\b",
                                            116: "\t",
                                            110: "\n",
                                            102: "\f",
                                            114: "\r"
                                        }, Y = function() {
                                            throw L = B = null,
                                                c()
                                        }, U = function() {
                                            for (var t, e, n, r, i, o = B, a = o.length; L < a; )
                                                switch (i = o.charCodeAt(L)) {
                                                    case 9:
                                                    case 10:
                                                    case 13:
                                                    case 32:
                                                        L++;
                                                        break;
                                                    case 123:
                                                    case 125:
                                                    case 91:
                                                    case 93:
                                                    case 58:
                                                    case 44:
                                                        return t = j ? o.charAt(L) : o[L],
                                                            L++,
                                                            t;
                                                    case 34:
                                                        for (t = "@",
                                                                L++; L < a; )
                                                            if (i = o.charCodeAt(L),
                                                            i < 32)
                                                                Y();
                                                            else if (92 == i)
                                                                switch (i = o.charCodeAt(++L)) {
                                                                    case 92:
                                                                    case 34:
                                                                    case 47:
                                                                    case 98:
                                                                    case 116:
                                                                    case 110:
                                                                    case 102:
                                                                    case 114:
                                                                        t += F[i],
                                                                            L++;
                                                                        break;
                                                                    case 117:
                                                                        for (e = ++L,
                                                                                n = L + 4; L < n; L++)
                                                                            i = o.charCodeAt(L),
                                                                            i >= 48 && i <= 57 || i >= 97 && i <= 102 || i >= 65 && i <= 70 || Y();
                                                                        t += M("0x" + o.slice(e, L));
                                                                        break;
                                                                    default:
                                                                        Y()
                                                                }
                                                            else {
                                                                if (34 == i)
                                                                    break;
                                                                for (i = o.charCodeAt(L),
                                                                        e = L; i >= 32 && 92 != i && 34 != i; )
                                                                    i = o.charCodeAt(++L);
                                                                t += o.slice(e, L)
                                                            }
                                                        if (34 == o.charCodeAt(L))
                                                            return L++,
                                                                t;
                                                        Y();
                                                    default:
                                                        if (e = L,
                                                        45 == i && (r = !0,
                                                            i = o.charCodeAt(++L)),
                                                        i >= 48 && i <= 57) {
                                                            for (48 == i && (i = o.charCodeAt(L + 1),
                                                            i >= 48 && i <= 57) && Y(),
                                                                    r = !1; L < a && (i = o.charCodeAt(L),
                                                            i >= 48 && i <= 57); L++)
                                                                ;
                                                            if (46 == o.charCodeAt(L)) {
                                                                for (n = ++L; n < a && (i = o.charCodeAt(n),
                                                                i >= 48 && i <= 57); n++)
                                                                    ;
                                                                n == L && Y(),
                                                                    L = n
                                                            }
                                                            if (i = o.charCodeAt(L),
                                                            101 == i || 69 == i) {
                                                                for (i = o.charCodeAt(++L),
                                                                    43 != i && 45 != i || L++,
                                                                        n = L; n < a && (i = o.charCodeAt(n),
                                                                i >= 48 && i <= 57); n++)
                                                                    ;
                                                                n == L && Y(),
                                                                    L = n
                                                            }
                                                            return +o.slice(e, L)
                                                        }
                                                        if (r && Y(),
                                                        "true" == o.slice(L, L + 4))
                                                            return L += 4,
                                                                !0;
                                                        if ("false" == o.slice(L, L + 5))
                                                            return L += 5,
                                                                !1;
                                                        if ("null" == o.slice(L, L + 4))
                                                            return L += 4,
                                                                null;
                                                        Y()
                                                }
                                            return "$"
                                        }, D = function(t) {
                                            var e, n;
                                            if ("$" == t && Y(),
                                            "string" == typeof t) {
                                                if ("@" == (j ? t.charAt(0) : t[0]))
                                                    return t.slice(1);
                                                if ("[" == t) {
                                                    for (e = []; t = U(),
                                                    "]" != t; n || (n = !0))
                                                        n && ("," == t ? (t = U(),
                                                        "]" == t && Y()) : Y()),
                                                        "," == t && Y(),
                                                            e.push(D(t));
                                                    return e
                                                }
                                                if ("{" == t) {
                                                    for (e = {}; t = U(),
                                                    "}" != t; n || (n = !0))
                                                        n && ("," == t ? (t = U(),
                                                        "}" == t && Y()) : Y()),
                                                        "," != t && "string" == typeof t && "@" == (j ? t.charAt(0) : t[0]) && ":" == U() || Y(),
                                                            e[t.slice(1)] = D(U());
                                                    return e
                                                }
                                                Y()
                                            }
                                            return t
                                        }, X = function(t, e, n) {
                                            var r = Z(t, e, n);
                                            r === y ? delete t[e] : t[e] = r
                                        }, Z = function(t, e, n) {
                                            var r, i = t[e];
                                            if ("object" == typeof i && i)
                                                if (g.call(i) == O)
                                                    for (r = i.length; r--; )
                                                        X(i, r, n);
                                                else
                                                    _(i, function(t) {
                                                        X(i, t, n)
                                                    });
                                            return n.call(t, e, i)
                                        };
                                        e.parse = function(t, e) {
                                            var n, r;
                                            return L = 0,
                                                B = "" + t,
                                                n = D(U()),
                                            "$" != U() && Y(),
                                                L = B = null,
                                                e && g.call(e) == w ? Z((r = {},
                                                    r[""] = n,
                                                    r), "", e) : n
                                        }
                                    }
                                }
                                return e.runInContext = o,
                                    e
                            }
                            var a = n(29)
                                , s = {
                                function: !0,
                                object: !0
                            }
                                , u = s[typeof e] && e && !e.nodeType && e
                                , l = s[typeof window] && window || this
                                , c = u && s[typeof t] && t && !t.nodeType && "object" == typeof i && i;
                            if (!c || c.global !== c && c.window !== c && c.self !== c || (l = c),
                            u && !a)
                                o(l, u);
                            else {
                                var f = l.JSON
                                    , h = l.JSON3
                                    , d = !1
                                    , p = o(l, l.JSON3 = {
                                    noConflict: function() {
                                        return d || (d = !0,
                                            l.JSON = f,
                                            l.JSON3 = h,
                                            f = h = null),
                                            p
                                    }
                                });
                                l.JSON = {
                                    parse: p.parse,
                                    stringify: p.stringify
                                }
                            }
                            a && (r = function() {
                                return p
                            }
                                .call(e, n, e, t),
                                !(void 0 !== r && (t.exports = r)))
                        }
                    ).call(this)
                }
            ).call(e, n(28)(t), function() {
                return this
            }())
        }
        , function(t, e) {
            t.exports = function(t) {
                return t.webpackPolyfill || (t.deprecate = function() {}
                    ,
                    t.paths = [],
                    t.children = [],
                    t.webpackPolyfill = 1),
                    t
            }
        }
        , function(t, e) {
            (function(e) {
                    t.exports = e
                }
            ).call(e, {})
        }
        , function(t, e, n) {
            "use strict";
            var r = n(31)
                , i = n(32)
                , o = {
                version: "0.2.2"
            }
                , a = function() {
                this.actionTypes = {},
                    this.storeQueue = [],
                    this.id = Date.now() + Math.round(1e3 * Math.random()),
                    this.middlewareQueue = new r(function(t) {
                        this.__invokeCallback__(t)
                    }
                        .bind(this),(!0))
            };
            a.prototype = {
                __invokeCallback__: function(t) {
                    this.storeQueue.forEach(function(e) {
                        var n, r = e.callbacks[t.type];
                        "function" == typeof r && (n = r(e.store, t),
                        void 0 !== n && e.store.event.publish(n))
                    })
                },
                use: function(t) {
                    "function" == typeof t && this.middlewareQueue.enter(t)
                },
                __dispatch__: function(t, e) {
                    var n, r = e(), i = this.actionTypes, o = r.type;
                    if (!o)
                        throw new Error("action type does not exist in \n" + JSON.stringify(r, null, 2));
                    if (n = i[o]) {
                        if (n !== t)
                            throw new Error('action type "' + o + '" is duplicate')
                    } else
                        i[o] = t;
                    this.middlewareQueue.execute(r)
                },
                createActions: function(t) {
                    var e, n, r = (this.id++).toString(32), i = this, o = {};
                    for (e in t)
                        n = t[e],
                            o[e] = function(t, e) {
                                return function() {
                                    var n = arguments;
                                    i.__dispatch__(e, function() {
                                        return t.apply(null, Array.prototype.slice.call(n))
                                    })
                                }
                            }(n, r);
                    return o
                },
                createMutableStore: function(t, e) {
                    if (!e)
                        throw new Error("schema must in createMutableStore arguments");
                    var n = new i(t)
                        , r = {
                        mutable: {},
                        event: {}
                    };
                    return r.mutable.get = n.mutable.get.bind(n.mutable),
                        r.event.subscribe = n.event.subscribe.bind(n.event),
                        r.event.unsubscribe = n.event.unsubscribe.bind(n.event),
                        this.storeQueue.push({
                            store: n,
                            callbacks: e
                        }),
                        r
                }
            },
                o.Dispatcher = a,
                t.exports = o
        }
        , function(t, e) {
            "use strict";
            var n = function(t, e) {
                this.workflows = [],
                    this.completeCallback = t,
                e && (this._workflows = [])
            };
            n.prototype = {
                enter: function(t) {
                    this.workflows.push(t),
                    this._workflows && this._workflows.push(t)
                },
                execute: function(t, e) {
                    e = e || this.workflows.concat();
                    var n;
                    e.length ? (n = e.shift())(t, this.execute.bind(this, t, e)) : (this._workflows && (this.workflows = this._workflows.concat()),
                        e = null,
                        this.completeCallback(t))
                }
            },
                t.exports = n
        }
        , function(t, e, n) {
            "use strict";
            var r = n(33)
                , i = n(34)
                , o = {
                string: !0,
                number: !0,
                null: !0,
                undefind: !0,
                boolean: !0
            }
                , a = function(t) {
                this.store = {},
                    Object.keys(t).forEach(function(e) {
                        this.store[e] = t[e]
                    }
                        .bind(this))
            };
            a.prototype = {
                set: function(t, e) {
                    if (t in this.store)
                        return this.store[t] = e,
                            t
                },
                get: function(t) {
                    var e = this.store[t]
                        , n = typeof e;
                    return o[n] ? e : r(e)
                },
                delete: function(t) {
                    return delete this.store[t],
                        t
                }
            };
            var s = function(t) {
                this.mutable = new a(t),
                    this.event = new i
            };
            t.exports = s
        }
        , function(t, e) {
            function n(t, e, r) {
                if (!(t instanceof Object))
                    return t;
                var i, o = Object.prototype.toString.call(t).slice(8, -1);
                switch (o) {
                    case "Array":
                        i = [];
                        break;
                    case "Date":
                        i = new Date(t.getTime());
                        break;
                    case "RegExp":
                        i = new RegExp(t);
                        break;
                    case "Function":
                        break;
                    case "Uint8Array":
                    case "Uint8ClampedArray":
                    case "Uint16Array":
                    case "Uint32Array":
                    case "Int8Array":
                    case "Int16Array":
                    case "Int32Array":
                    case "Float32Array":
                    case "Float64Array":
                        i = t.subarray();
                        break;
                    default:
                        i = {}
                }
                if (e.push(t),
                    r.push(i),
                t instanceof Array)
                    for (var a = 0; a < t.length; a++)
                        i[a] = n(t[a], e, r);
                for (var s = Object.keys(t).sort(), u = Object.keys(i).sort(), l = 0; l < s.length; l++) {
                    var c = s[l];
                    if (u.length > 0 && c === u[0])
                        u.shift();
                    else if (Object.prototype.hasOwnProperty.call(t, c)) {
                        var f = t[c]
                            , h = e.indexOf(f);
                        i[c] = h !== -1 ? r[h] : n(t[c], e, r)
                    }
                }
                return i
            }
            t.exports = function(t) {
                return n(t, [], [])
            }
        }
        , function(t, e) {
            "use strict";
            var n = function() {
                this.handlers = []
            };
            n.prototype = {
                publish: function(t) {
                    this.handlers.forEach(function(e) {
                        e.type ? e.type === t && e.handler(t) : e.handler(t)
                    })
                },
                subscribe: function(t, e) {
                    var n = {};
                    "function" == typeof t ? n.handler = t : (n.handler = e,
                        n.type = t);
                    for (var r, i = 0; i < this.handlers.length; i++)
                        r = this.handlers[i],
                        r.type === t && this.handlers.splice(i, 1);
                    this.handlers.push(n)
                },
                unsubscribe: function(t, e) {
                    "function" == typeof t && (e = t,
                        t = null);
                    for (var n, r = 0, i = !1; r < this.handlers.length; r++)
                        n = this.handlers[r],
                            n.type ? t && e ? i = n.type === t && n.handler === e : t ? i = n.type === t : e && (i = n.handler === e) : i = n.handler === e,
                        i && this.handlers.splice(r--, 1)
                }
            },
                t.exports = n
        }
        , function(t, e) {
            "use strict";
            function n() {
                window.yoda_doc.innerHTML = "",
                    window.seed[window.yoda_callModule].initModule[window.yoda_callMethod]()
            }
            Object.defineProperty(e, "__esModule", {
                value: !0
            });
            var r = "<div style='position: fixed; background: rgba(0, 0, 0, 0.4); width: 100%; height: 100%; left: 0; right: 0; top: 0; bottom: 0; z-index: 997;'>\n        <div style='position: fixed; left: 3em; right: 3em; bottom: 50%; border-radius: 10px; background: #FCFCFC; z-index: 998; height:auto;' id='childRoot'></div>\n    </div>"
                , i = "<div style='position: fixed; background: rgba(0, 0, 0, 0.4); width: 100%; height: 100%; left: 0; right: 0; top: 0; bottom: 0; z-index: 997;'>\n        <div style='position: fixed; left: 3em; right: 3em; height: 100%; border-radius: 10px; background: rgba(0, 0, 0, 0); z-index: 998;' id='childRoot'></div>\n    </div>"
                , o = window.seed.env || "pro"
                , a = function(t) {
                var e = t.root;
                window.yoda_doc = e,
                    e.innerHTML = r,
                    t.root = "childRoot",
                    window.yoda_callModule = t.succModule,
                    window.yoda_callMethod = t.succCallbackFun,
                    window.yoda_executeGlobalCallBack = n,
                    t.succCallbackFun = "yoda_executeGlobalCallBack",
                    window.YodaSeed(t, o)
            }
                , s = function(t) {
                var e = t.root;
                window.yoda_doc = e,
                    e.innerHTML = i,
                    t.root = "childRoot",
                    window.yoda_callModule = t.succModule,
                    window.yoda_callMethod = t.succCallbackFun,
                    window.yoda_executeGlobalCallBack = n,
                    t.succCallbackFun = "yoda_executeGlobalCallBack",
                    window.YodaSeed(t, o)
            }
                , u = {
                iAdapter: a,
                pcAdapter: s
            };
            e.default = u
        }
    ]);
    return {
        babelHelpers: babelHelpers,
        window: window,
        navigator: navigator,
        document: document,
        screen: screen
    }
}

exports.abcdefg = abcdefg