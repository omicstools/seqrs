# GitHub Actions 发布流程说明

## 自动CI/CD流水线

该项目已配置GitHub Actions自动化流程：

### 1. CI 流程 (ci.yml)
在以下情况自动触发：
- 推送到 `main`, `master`, `develop` 分支
- 创建Pull Request到这些分支

**执行步骤：**
- 运行测试 (`cargo test`)
- 构建发布二进制 (`cargo build --release`)
- 上传Linux amd64二进制到Artifacts

### 2. Release 流程 (release.yml)
在创建版本标签时自动触发：
- 创建形如 `v*` 的标签（例如 `v1.0.0`, `v1.1.0`）

**执行步骤：**
- 构建发布二进制
- 创建压缩包及SHA256校验和
- 自动创建GitHub Release并上传文件

## 如何发布新版本

### 方式一：使用Git命令（推荐）

```bash
# 1. 更新版本号在 Cargo.toml
# 修改 version = "1.x.x"

# 2. 更新 CHANGELOG.md
# 在 [Unreleased] 部分添加变更内容，然后创建新的 [x.y.z] 段落

# 3. 提交更改
git add Cargo.toml CHANGELOG.md
git commit -m "chore: bump version to 1.x.x"

# 4. 创建标签并推送
git tag -a v1.x.x -m "Release version 1.x.x"
git push origin main
git push origin v1.x.x
```

### 方式二：使用GitHub Web界面

1. 在GitHub仓库页面，点击 "Releases"
2. 点击 "Draft a new release"
3. 创建新标签（格式：v1.x.x）
4. 填写Release标题和说明
5. 点击 "Publish release"
6. GitHub Actions 会自动构建和上传文件

## 版本号规范

遵循 [语义化版本控制](https://semver.org/lang/zh-CN/)：

- **主版本号 (Major)**: 不兼容的API变更
- **次版本号 (Minor)**: 向下兼容的功能性新增
- **修订号 (Patch)**: 向下兼容的问题修正

例如：
- `v1.0.0` - 第一个发布版本
- `v1.1.0` - 新增功能
- `v1.0.1` - Bug修复
- `v2.0.0` - 主要版本升级

## 发布包内容

自动生成的发布包包含：
- `seqrs-vX.Y.Z-linux-amd64.tar.gz` - 压缩的二进制和文档
- `seqrs-vX.Y.Z-linux-amd64.tar.gz.sha256` - SHA256校验和文件

## 校验下载的二进制

```bash
# 解压
tar -xzf seqrs-v1.0.0-linux-amd64.tar.gz

# 验证校验和
sha256sum -c seqrs-v1.0.0-linux-amd64.tar.gz.sha256

# 运行程序
./seqrs-v1.0.0-linux-amd64/seqrs --help
```

## 故障排除

### Release流程失败
- 检查GitHub Actions日志：仓库 → Actions → 找到失败的工作流
- 确保标签格式正确（`v` 开头，例如 `v1.0.0`）
- 检查Cargo.toml语法是否正确

### 二进制无法执行
- 确认系统架构为amd64 (x86-64)
- 检查文件权限：`chmod +x seqrs`
- 确认依赖库已安装（Linux glibc 2.29+）
