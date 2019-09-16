let log4js = require('log4js')
let path = require('path')

let programName = 'data_cover'
let dirName = new Date().toLocaleDateString().toString().replace(/\//g, '-')
let filePath = path.join(__dirname, 'log/dataCover', dirName, programName)

log4js.configure({
  appenders: {
    // 记录器1:输出到控制台
    console: {
      type: 'console'
    },
    // 记录器3：输出到日期文件
    dateFile: {
      type: 'dateFile',
      // 您要写入日志文件的路径
      filename: `${filePath}`,
      // 时间文件 保存多少天，距离当前天daysToKeep以前的log将被删除
      daysToKeep: 30,
      // maxLogSize: 10,
      alwaysIncludePattern: true,
      // compress : true,//（默认为false） - 在滚动期间压缩备份文件（备份文件将具有.gz扩展名）
      // （可选，默认为.yyyy-MM-dd） - 用于确定何时滚动日志的模式。格式:.yyyy-MM-dd-hh:mm:ss.log
      pattern: '_hh.log',
      // default "utf-8"，文件的编码
      encoding: 'utf-8'
    },
    sqlLog: {
      type: 'dateFile',
      // 您要写入日志文件的路径
      filename: `${filePath.replace(/data_cover/, 'data_sql')}`,
      // （默认为false） - 将模式包含在当前日志文件的名称以及备份中
      alwaysIncludePattern: true,
      // （可选，默认为.yyyy-MM-dd） - 用于确定何时滚动日志的模式。格式:.yyyy-MM-dd-hh:mm:ss.log
      pattern: '_hh.log',
      daysToKeep: 30,
      // maxLogSize: 10,
      // default "utf-8"，文件的编码
      encoding: 'utf-8'
    }
  },
  categories: {
    // 默认log类型，输出到控制台 log文件 log日期文件 且登记大于info即可
    default: { appenders: ['dateFile', 'console'], level: 'info' },
    sql_log: { appenders: ['sqlLog'], level: 'info' }
  }
})

let s = log4js.getLogger('sql_log')
s.info(`
SELECT
count( * ) soldNum 
FROM
wxapp.shop_order_2 
WHERE
( originId = '4937020' OR shopId = '4937020' ) 
AND completeTime > '1559750399' 
AND completeTime < '1562342399';
`)

// module.exports = log4js.getLogger('console')