var CryptoJS = require("crypto-js");

module.exports = class {
    constructor(key, iv) {
        this.key = key
        this.iv = iv
    }

    // 加密
    encrypt_aes_pkcs(data) {
        let key = CryptoJS.enc.Utf8.parse(this.key);
        let iv = CryptoJS.enc.Utf8.parse(this.iv);

        var encrypted = CryptoJS.AES.encrypt(data, key, {
            iv: '',
            mode: CryptoJS.mode.ECB,
            padding: CryptoJS.pad.Pkcs7
        });
        // 转换为字符串
        encrypted = encrypted.toString();

        return encrypted
    }

    // 解密
    decrypt_aes_pkcs(data) {
        var decrypted = CryptoJS.AES.decrypt(data, key, {
            iv: '',
            mode: CryptoJS.mode.ECB,
            padding: CryptoJS.pad.Pkcs7
        });
    
        // 转换为 utf8 字符串
        decrypted = CryptoJS.enc.Utf8.stringify(decrypted);
        return decrypted
    }
}

// let key = "jvzempodf8f9anyt";
// let iv = "jvzempodf8f9anyt";
// let _aes = new AESCipher(key, iv)

// var x = JSON.stringify({"ts":1553954578591,"cts":1553954719379,"brVD":[375,667],"brR":[[375,667],[375,667],24,24],"aM":""})
// var a = _aes.encrypt_aes_pkcs(x)
// console.log(a)
// console.log(_aes.decrypt_aes_pkcs(a))