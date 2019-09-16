const router = require('koa-router')()
const AESCipher = require('../layout/v2/x_for_with')


// 美团 X-FOR-WITH 加密
router.post('/aes_pkcs', async (ctx, next) => {
    const { request } = ctx;
    let data = request.body;
    let _aes = new AESCipher(data.key, data.iv)
    let encrypt = _aes.encrypt_aes_pkcs(data.data)
    ctx.response.body = encrypt
})

module.exports = router