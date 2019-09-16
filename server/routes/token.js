const router = require('koa-router')();
const ht_slider = require('../layout/v2/ht_slider');
const Rohr_Opt = require('../layout/v2/token')
const m_Rohr_Opt = require('../layout/v2/m_token')
const m_info_slider = require('../layout/v2/m_info_slider');
const dp_token = require('../layout/v2/dp_slider');

router.get('/', async (ctx, next) => {
  ctx.response.body = 'hello token!'
});

// 大众点评滑块加密
router.post('/dpToken', async (ctx, next) => {
  const config = ctx.request.body;

  config.country = '中国大陆';
  config.defaultIndex = Number(config.defaultIndex);

  for (let key in config) {
    if (config[key] === 'null') {
      config[key] = null
    } else if (config[key] === 'false' || config[key] === 'False') {
      config[key] = false
    }
  }

  let dpToken = dp_token.get_behavior_token(config);
  ctx.response.body = dpToken
})

// 美团 m 端 详情页 滑块参数加密逻辑
router.post('/m_info_token', async (ctx, next) => {
  const data = ctx.request.body;

  let location = typeof data.location == 'string' ? JSON.parse(data.location) : data.location;
  let user_agent = data.user_agent;
  let config = data.config

  config.country = '中国大陆';
  config.defaultIndex = Number(config.defaultIndex);

  for (let key in config) {
    if (config[key] === 'null') {
      config[key] = null
    } else if (config[key] === 'false' || config[key] === 'False') {
      config[key] = false
    }
  }

  let m_info_token = m_info_slider.get_behavior_token(config, location, String(user_agent))
  ctx.response.body = m_info_token
})

router.post('/m_token', async (ctx, next) => {
  const { request } = ctx;
  let data = request.body;
  let m_token = m_Rohr_Opt.Rohr_Opt.reload(data)
  ctx.response.body = m_token
});

// 美团商家后台登陆token加密逻辑
router.post('/token', async (ctx, next) => {
  // "https://epassport.meituan.com/api/account/login?loginContinue=http://e.waimai.meituan.com/v2/epassport/entry&&only_auth=undefined"
  let token = Rohr_Opt.Rohr_Opt.reload("https://epassport.meituan.com/api/account/login?loginContinue=http://e.waimai.meituan.com/v2/epassport/entry&&only_auth=undefined")
  ctx.response.body = token
});

// 美团商家后台登陆滑块加密逻辑
router.post('/tbToken', async (ctx, next) => {
  const { request } = ctx;
  let config = request.body;
  
  config.country = '中国大陆';
  config.defaultIndex = Number(config.defaultIndex);

  for (let key in config) {
    if (config[key] === 'null') {
      config[key] = null
    } else if (config[key] === 'false' || config[key] === 'False') {
      config[key] = false
    }
  }

  let params = ht_slider.get_behavior_token(config);

  ctx.response.body = params
})

// 美团商家后台登陆图片验证token加密逻辑
router.post('/imgToken', async (ctx, next) => {
  const { request } = ctx;

  let { request_code, result_code } = request.body

  let params = ht_slider.get_img_token(request_code, result_code)

  ctx.response.body = params
})

module.exports = router
